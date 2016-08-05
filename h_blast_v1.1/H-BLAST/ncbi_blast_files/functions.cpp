#define CURRENT(a) ((a)&0x1)
#define NEXT(a) CURRENT(a+1)
#define MAX_CPU_CORE_NO 100
#define MIN(a,b) ((a>b)?b:a)
#define INT_CEIL(a,b) ((a+b-1)/b) /// The ceiling function of an integer r.w.t. an integer b.

#include <math.h>


static double max_cpu_walltime[3];
static double current_cpu_walltime[MAX_CPU_CORE_NO];

static int load_balancer_flag = 1, load_factor=4;

static double gpu_cpu_capability_rate=12;

static float Gpu_db_used_ratio_in_byte=-1;

static char mem_ecc_flag = 0;

static char flag_fixed_gpu_db_used_ratio=1;

static int subject_chunck_size = 200;

static int base_value = 100;
static unsigned int length_per_query = 1000;


void init_db_volume_size_in_byte(unsigned long* db_volume_size_in_byte)
{
    if(NULL==db_volume_size_in_byte)
    {
        printf("Error: the pointer in \"init_db_volume_size_in_byte\" is NULL!\n ");
        gpu_abnormally_exit(-1);
    }
    int i;
    for(i=0; i<3; i++)
        db_volume_size_in_byte[i] = 0;
}

void init_max_and_current_cpu_walltime()
{
    int i;
    for(i=0; i<MAX_CPU_CORE_NO; i++)
        current_cpu_walltime[i]=0;
    for(i=0; i<3; i++)
        max_cpu_walltime[i] = 0;
}

float get_new_gpudb_ratio(const unsigned int volume_id, const float old_gpudb_ratio, const unsigned int no_of_cpu_threads, const unsigned int no_of_gpu_cards,
                          const unsigned long* db_volume_size_in_byte, const float* old_gpudb_ratio_array, const long int * old_successful_extension_array,
                          const Boolean is_gapped_calculation, const int query_length)
{
    double gpu_walltime = get_gpu_walltime();
    double cpu_walltime = max_cpu_walltime[((volume_id-1)%3)];

    double cpu_db_size_ratio, gpu_db_size_ratio;
    if( (0==db_volume_size_in_byte[0]) || (0==db_volume_size_in_byte[1]) ||(0==db_volume_size_in_byte[2]) )
    {
        printf("The db size in byte is 0!\n");
        gpu_abnormally_exit(-1);
    }
    if(0.0==old_gpudb_ratio_array[0] || 0.0==old_gpudb_ratio_array[1])
    {
        printf("The rate for gpu is 0!\n");
        gpu_abnormally_exit(-1);
    }


    double rate1 = (double)old_gpudb_ratio_array[CURRENT(volume_id)]*0.96;
    double rate2 = (double)old_gpudb_ratio_array[CURRENT(volume_id-1)]*0.96;



    if(is_gapped_calculation != 0)
    {

        rate1 = (double)(db_volume_size_in_byte[(volume_id)%3]) * (1.0 - rate1) + (double)old_successful_extension_array[CURRENT(volume_id)]*(base_value + (double)(750-base_value)*MIN(1500, length_per_query)/1500 )/(rate1+0.03) ;
        rate2 = (double)(db_volume_size_in_byte[(volume_id-1)%3]) * (1.0-rate2) + (double)old_successful_extension_array[CURRENT(volume_id-1)]*(base_value + (double)(750-base_value)*MIN(1500,length_per_query)/1500 )/(rate2+0.03);
    }
    else
    {

        if(load_balancer_flag==0)
        {
            rate1 = (double)(db_volume_size_in_byte[(volume_id)%3]) * (1.0 - rate1);
            rate2 = (double)(db_volume_size_in_byte[(volume_id-1)%3]) * (1.0-rate2);
        }
        else
        {
            rate1 = (double)(db_volume_size_in_byte[(volume_id)%3]) * (1.0 - rate1) + (double)old_successful_extension_array[CURRENT(volume_id)]*load_factor/(rate1+0.03) ;
            rate2 = (double)(db_volume_size_in_byte[(volume_id-1)%3]) * (1.0-rate2) + (double)old_successful_extension_array[CURRENT(volume_id-1)]*load_factor/(rate1+0.03) ;
        }

    }


    cpu_db_size_ratio = rate1 / rate2;
    gpu_db_size_ratio = (double)(db_volume_size_in_byte[(volume_id+1)%3]) / db_volume_size_in_byte[(volume_id)%3];


    if( (query_length>1000) && (query_length<4000) )
        gpu_db_size_ratio = pow(gpu_db_size_ratio, 0.9);
    else if( (query_length>=4000) && (query_length<8000) )
        gpu_db_size_ratio = pow(gpu_db_size_ratio, 0.8);
    else if(query_length>=8000)
        gpu_db_size_ratio = pow(gpu_db_size_ratio, 0.75);


    cpu_walltime *= cpu_db_size_ratio;
    gpu_walltime *= gpu_db_size_ratio;


    float new_rate=0;

    if(0==gpu_walltime) return 0.9;

    float equivalent_no_of_gpu_cards = no_of_cpu_threads / gpu_cpu_capability_rate + no_of_gpu_cards;

    float adjust_ratio = (cpu_walltime - gpu_walltime) / gpu_walltime;


    const double E_RATIO_OF_DELTA_TIME=0.65;
    double M_rate = (double)(E_RATIO_OF_DELTA_TIME*no_of_gpu_cards);
    M_rate = M_rate/(M_rate+(1-E_RATIO_OF_DELTA_TIME)*no_of_cpu_threads / gpu_cpu_capability_rate);

    if((cpu_walltime - gpu_walltime)>0)
        adjust_ratio *= M_rate/equivalent_no_of_gpu_cards;
    else adjust_ratio *= 1.0 - M_rate/equivalent_no_of_gpu_cards;


    new_rate = old_gpudb_ratio + adjust_ratio/0.96;

    new_rate = (new_rate>1)?1:new_rate;

    new_rate = (new_rate<0)?0:new_rate;

    return new_rate;
}

void Read_H_BLAST_runtime_options()
{
    FILE * pFile;

    char* gpu_db_used_ratio_file_name = "./H-BLAST_runtime_options";

    pFile = fopen ( gpu_db_used_ratio_file_name, "rb" );

    if (pFile==NULL)
    {
        printf ("Warning: There is no user defined Gpu db used ratio! Use the default one!\n");
        return;
    }


    char *line_read = (char *) malloc( sizeof(char) );
    size_t nbytes = 1;
    int bytes_read, max_storage_tmp = 0;

    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    Gpu_db_used_ratio_in_byte = ( (atof(line_read) < 0 ) || (atof(line_read) >1 ) ) ? 0.0 :(double)atof(line_read);

    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    flag_fixed_gpu_db_used_ratio = ( 1 == atoi(line_read) ) ? 1:0;

    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    gpu_cpu_capability_rate = ( atof(line_read) <= 0 ) ? 12 :(double)atof(line_read);

    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    max_storage_tmp = ( atoi(line_read) > 0 ) ? atoi(line_read) : 20 ;
    set_MAX_SUCCESSFUL_HITS_PER_SEQUENCE(max_storage_tmp);

    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    mem_ecc_flag = ( atoi(line_read) != 0 ) ? 1 : 0 ;


    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    load_balancer_flag = ( atoi(line_read) != 0 ) ? 1 : 0 ;

    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    load_factor = ( atoi(line_read) < 0 ) ? 4 :  atoi(line_read) ;

    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    subject_chunck_size = ( atoi(line_read) < 0 ) ? 200 :  atoi(line_read) ;

    bytes_read = getline(&line_read, &nbytes, pFile);
    bytes_read = getline(&line_read, &nbytes, pFile);
    base_value = ( atoi(line_read) < 0 ) ? 100 :  atoi(line_read) ;

    free(line_read);
    fclose(pFile);
}


void get_GPU_result( Int4** _h_Hits,
                     const Int4 _h_Hits_bytes,
                     GPUBlastInitHitList** _h_GPUBlastInitHitList,
                     const Int4 _h_GPUBlastInitHitList_bytes,
                     const Int4 _num_sequences,
                     double * _elapsed_time)
{
    Int4* h_Hits_temp=NULL;
    GPUBlastInitHitList* h_GPUBlastInitHitList_temp=NULL;

    struct timespec CPU_start, CPU_end;


    h_Hits_temp = (Int4*) realloc(*_h_Hits, _h_Hits_bytes);

    if( NULL != h_Hits_temp )
        *_h_Hits = h_Hits_temp;
    else
    {
        puts("Error (re)allocating memory for the GPU\n");
        all_gpu_abnormally_exit(1);
    }

    h_GPUBlastInitHitList_temp = (GPUBlastInitHitList*) realloc(*_h_GPUBlastInitHitList, _h_GPUBlastInitHitList_bytes);
    if( NULL != h_GPUBlastInitHitList_temp )
        *_h_GPUBlastInitHitList = h_GPUBlastInitHitList_temp;
    else
    {
        printf("Error (re)allocating memory for the GPU\n");
        all_gpu_abnormally_exit(1);
    }

    clock_gettime(CLOCK_MONOTONIC, &CPU_start);

    H_BLAST_get_data(*_h_Hits,  _h_Hits_bytes,
                        *_h_GPUBlastInitHitList,
                        _h_GPUBlastInitHitList_bytes,
                        _num_sequences);

    clock_gettime(CLOCK_MONOTONIC, &CPU_end);
    if(_elapsed_time ==NULL)
    {
        fprintf(stderr, "get_GPU_result: Error, \"_elapsed_time\" is invaild!\n");
        all_gpu_abnormally_exit(1);
    }
    *_elapsed_time = (CPU_end.tv_sec - CPU_start.tv_sec) + (double)(CPU_end.tv_nsec - CPU_start.tv_nsec) / BILLION;

    return;

}


void Write_GPU_information_file_preliminary(FILE* GPU_information_file,
        const Int4 num_blocksx,
        const Int4 num_threadsx,
        const char percentage,
        const Int4 num_volumes)
{
    fprintf(GPU_information_file,"#Number of blocksx\n");
    fprintf(GPU_information_file,"%d\n", num_blocksx );
    fprintf(GPU_information_file,"#Number of threadsx\n");
    fprintf(GPU_information_file,"%d\n", num_threadsx);
    fprintf(GPU_information_file,"#Percentage\n");
    fprintf(GPU_information_file,"%d\n", percentage);
    fprintf(GPU_information_file,"#Number of volumes\n");
    fprintf(GPU_information_file,"%d\n", num_volumes );

}



void Write_GPU_information_file(FILE* GPU_information_file,
                                const char* GPU_Database_name,
                                const unsigned long H_BLAST_database_bytes,
                                const Int4 volume_length,
                                const Int4 num_sequences_to,
                                const Int4 store_limit,
                                const Int4 break_limit,
                                const Int4 Group_number,
                                const Int4* Sequence_length_vector,
                                const unsigned long CPU_rest_bytes)
{

    fprintf(GPU_information_file,"#Volume name\n");
    fprintf(GPU_information_file,"%s\n", GPU_Database_name);
    fprintf(GPU_information_file,"#Volume size in bytes\n");
    fprintf(GPU_information_file,"%ld\n", H_BLAST_database_bytes );
    fprintf(GPU_information_file,"#Volume sequences\n");
    fprintf(GPU_information_file,"%d\n", volume_length );
    fprintf(GPU_information_file,"#Number of sequences in GPU volume\n");
    fprintf(GPU_information_file,"%d\n", num_sequences_to );
    fprintf(GPU_information_file,"#Store limit\n");
    fprintf(GPU_information_file,"%d\n", store_limit );
    fprintf(GPU_information_file,"#Break limit\n");
    fprintf(GPU_information_file,"%d\n", break_limit );
    fprintf(GPU_information_file,"#Group number\n");
    fprintf(GPU_information_file,"%d\n", Group_number );
    fprintf(GPU_information_file,"#Length of each group\n");
    Int4 i = 0;
    for(i = 0; i < Group_number; ++i)
        fprintf(GPU_information_file,"%d\n", Sequence_length_vector[i] );
    fprintf(GPU_information_file,"#Length of the rest sequences in CPU database\n");
    fprintf(GPU_information_file,"%ld\n",CPU_rest_bytes);

}


void Read_GPU_DB_information_file_preliminary(FILE** GPU_information_file,
        const BlastSeqSrc* seq_src,
        Int4* num_blocksx, Int4* num_threadsx,
        Int4* percentage, Int4* num_volumes)
{
    char* Database_name = BlastSeqSrcGetName(seq_src);
    char* GPU_information_name = (char*) calloc(strlen(Database_name) + strlen(".hDbI") + 1,sizeof(char));
    char* GPU_information_name1 = (char*) calloc(strlen(Database_name) + strlen(".gpuinfox") + 1,sizeof(char));
    strcat(GPU_information_name, Database_name);
    strcat(GPU_information_name, ".hDbI");
    strcat(GPU_information_name1, Database_name);
    strcat(GPU_information_name1, ".gpuinfox");

    *GPU_information_file = fopen( GPU_information_name, "rb");
    if( NULL == (*GPU_information_file) )
    {
        if( NULL == (*GPU_information_file = fopen( GPU_information_name1, "rb")) )
        {
        printf("GPU information file cannot be opened for reading. Exiting...\n");
        exit(1);
        }

    }

    char *line_read = (char *) malloc( sizeof(char) );
    size_t nbytes = 1;
    Int4 bytes_read = getline(&line_read, &nbytes, *GPU_information_file);
    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);
    *num_blocksx = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);
    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);
    *num_threadsx = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);
    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);
    *percentage = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);
    bytes_read = getline(&line_read, &nbytes, *GPU_information_file);
    *num_volumes = atoi(line_read);

    free(line_read);
    free(GPU_information_name);

}

void Read_GPU_DB_information_file(const FILE* GPU_information_file,
                               char** GPU_volume_name,
                               unsigned long* H_BLAST_database_bytes,
                               Int4* volume_length,
                               Int4* num_sequences_to,
                               Int4* store_limit,
                               Int4* break_limit,
                               Int4* Group_number,
                               unsigned long* db_size_in_byte,
                               Int4** Sequence_length_vector
                              )
{

    char *line_read = (char *) malloc( sizeof(char) );
    size_t nbytes = 1;
    Int4 bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read first line which is a comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the volume name
    char* GPU_volume_name_temp = (char*)realloc(*GPU_volume_name, bytes_read*sizeof(char) );
    if( NULL != GPU_volume_name_temp )
    {
        *GPU_volume_name = GPU_volume_name_temp;
        memset(*GPU_volume_name, 0, bytes_read*sizeof(char));
        memcpy(*GPU_volume_name, line_read, bytes_read - 1);
    }
    else
    {
        printf("Error (re)allocating memory for the GPU information file\n");
        exit(1);
    }

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the size of the volume in bytes
    *H_BLAST_database_bytes = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the number of sequences in this volume
    *volume_length = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the number of sequences processed by the GPU in this volume
    *num_sequences_to = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the ID of the last sequences processed by the GPU in this volume
    *store_limit = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the ID of the last sequences of this volume
    *break_limit = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the number of groups in this volume
    *Group_number = atoi(line_read);

    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    Int4* Sequence_length_vector_temp = (Int4*) realloc(*Sequence_length_vector, (*Group_number)*sizeof(Int4));
    if( NULL != Sequence_length_vector_temp )
        *Sequence_length_vector = Sequence_length_vector_temp;
    else
    {
        printf("Error (re)allocating memory for the GPU information file\n");
        exit(1);
    }
    Int4 i = 0;
    for( i = 0; i < *Group_number; ++i)
    {
        bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the length of each group of this volume
        (*Sequence_length_vector)[i] = atoi(line_read);
    }


    //Read CPU_rest_bytes
    unsigned long CPU_rest_bytes = 0;
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read comment
    bytes_read = getline(&line_read, &nbytes, GPU_information_file);//Read the number of groups in this volume
    CPU_rest_bytes = atoi(line_read);

    *db_size_in_byte = CPU_rest_bytes + *H_BLAST_database_bytes;

    free(line_read);

}


void get_GPU_DB_information(unsigned long* H_BLAST_database_bytes,
                            Int4* volume_length,
                            Int4* num_sequences_to,
                            Int4* store_limit,
                            Int4* break_limit,
                            Int4* Group_number,
                            Int4** Sequence_length_vector,
                            const Int4 cpu_threads,
                            const int thread_count_per_block,
                            const int gpu_db_group_size,
                            float* gpu_db_ratio_ptr )
{

    if(NULL==gpu_db_ratio_ptr)
    {
        printf("Error: gpu ratio is unaccessable!\n");
        exit(-1);
    }

    /// determind the size of GPU DB to be processed
    double ratio = 0;
    if(-1==(*gpu_db_ratio_ptr) )
    {
        switch ( cpu_threads )
        {
        case 1:
            ratio = 1;
            break;
        case 2:
            ratio = 0.9;
            break;
        case 3:
            ratio = 0.9;
            flag_fixed_gpu_db_used_ratio=0;
            break;
        case 4:
            ratio = 0.8;
            flag_fixed_gpu_db_used_ratio=0;
            break;
        case 5:
            ratio = 0.75;
            flag_fixed_gpu_db_used_ratio=0;
            break;
        case 6:
            ratio = 0.65;
            flag_fixed_gpu_db_used_ratio=0;
            break;
        case 7:
            ratio = 0.55;
            flag_fixed_gpu_db_used_ratio=0;
            break;
        case 8:
            ratio = 0.45;
            flag_fixed_gpu_db_used_ratio=0;
            break;
        default:
            ratio = 0.4;
            flag_fixed_gpu_db_used_ratio=0;
            break;

        }
        printf("new ratio is %lg\n", ratio);
    }
    else ratio = (*gpu_db_ratio_ptr);


    Int4 l=0;
    double total_workload_cost=0, actual_workload_cost=0;

    for(l=0; l<(*Group_number-1); l++)
        total_workload_cost += (*Sequence_length_vector)[l];
    total_workload_cost*=gpu_db_group_size;

    total_workload_cost += ((double)(*Sequence_length_vector)[*Group_number-1]) * (*num_sequences_to - (*Group_number-1)*gpu_db_group_size);
    actual_workload_cost = total_workload_cost * ratio;

    Int4 k = (*Group_number);
    double tmp_workload_cost = total_workload_cost;
    if(ratio<1)
    {
        tmp_workload_cost -= ((double)(*Sequence_length_vector)[*Group_number-1]) * (*num_sequences_to - (*Group_number-1)*gpu_db_group_size);
        k--;
        do
        {
            if(tmp_workload_cost<= actual_workload_cost )
                break;
            else
            {
                tmp_workload_cost -=  (double)(*Sequence_length_vector)[k-1]*gpu_db_group_size;
                k--;
            }

        }
        while(k>0);
    }


    if(k!=*Group_number)
    {
        double rest_workload_from_expected = actual_workload_cost - tmp_workload_cost;
        int no_threads_added = (0!=(*Sequence_length_vector)[k])? (rest_workload_from_expected / ( (*Sequence_length_vector)[k] *thread_count_per_block) )*thread_count_per_block : 0;
        tmp_workload_cost += (double)(*Sequence_length_vector)[k] *no_threads_added;
        if(flag_fixed_gpu_db_used_ratio!=1)
        {
            *gpu_db_ratio_ptr = tmp_workload_cost / total_workload_cost;
            *gpu_db_ratio_ptr = ((*gpu_db_ratio_ptr)>1)?1:(*gpu_db_ratio_ptr);
            *gpu_db_ratio_ptr = ((*gpu_db_ratio_ptr)<0)?0:(*gpu_db_ratio_ptr);
        }

        if( (k* gpu_db_group_size + no_threads_added) < *num_sequences_to)
        {
            *num_sequences_to = k* gpu_db_group_size + no_threads_added;
            *Group_number = k;
            *Group_number+= (0==no_threads_added)?0:1;

            unsigned long GPU_db_bytes_count = 0;
            for(l=0; l<k; l++)
                GPU_db_bytes_count += (*Sequence_length_vector)[l]*gpu_db_group_size;

            *H_BLAST_database_bytes = GPU_db_bytes_count + ((0!=no_threads_added)? (*Sequence_length_vector)[k]*INT_CEIL(no_threads_added, thread_count_per_block)*thread_count_per_block:0);
        }
    }
}

int ExtractVolumeSequences(unsigned char * buffer)
{
    //Extracting information from the index file (see index_file.txt for the format information)

    int version =
        buffer[0 + 0]*(1<<24) +
        buffer[0 + 1]*(1<<16) +
        buffer[0 + 2]*(1<<8) +
        buffer[0 + 3];
    int type =
        buffer[4 + 0]*(1<<24) +
        buffer[4 + 1]*(1<<16) +
        buffer[4 + 2]*(1<<8) +
        buffer[4 + 3];
    if( type != 1 )
    {
        fprintf(stderr,"Error index file is not a protein index file\n");
        exit(1);
    }
    int title_bytes =
        buffer[8 + 0]*(1<<24) +
        buffer[8 + 1]*(1<<16) +
        buffer[8 + 2]*(1<<8) +
        buffer[8 + 3];
    int date_bytes =
        buffer[8 + 4 + title_bytes + 0 ]*(1<<24) +
        buffer[8 + 4 + title_bytes + 1 ]*(1<<16) +
        buffer[8 + 4 + title_bytes + 2 ]*(1<<8) +
        buffer[8 + 4 + title_bytes + 3 ];
    int sequences =
        buffer[8 + 4 + 4 + title_bytes + date_bytes + 0 ]*(1<<24) +
        buffer[8 + 4 + 4 + title_bytes + date_bytes + 1 ]*(1<<16) +
        buffer[8 + 4 + 4 + title_bytes + date_bytes + 2 ]*(1<<8) +
        buffer[8 + 4 + 4 + title_bytes + date_bytes + 3 ];

    return sequences;
}

void ReadH_BLAST_volinfo(const BlastSeqSrc* seq_src, Int4 *num_volumes, Int4** volume_lengths)
{
    if( GPU_VERBOSE )
        printf("Reading volume information\n");

    Int4* volume_lengths_temp;
    *volume_lengths = NULL;

    //deep copy the database path to a temporary string because "dirname()" might modify the database path
    char* Database_name = BlastSeqSrcGetName(seq_src);
    int Database_name_bytes = strlen( Database_name ) + 1;
    char* Database_copy = (char*) calloc( Database_name_bytes, sizeof(char) );
    memcpy( Database_copy, Database_name, Database_name_bytes );

    char *aliasfile_name = (char*)calloc(strlen(Database_copy) + strlen(".pal") + 1, sizeof(char));
    strcat(aliasfile_name, Database_copy);
    strcat(aliasfile_name, ".pal");
    FILE* aliasfile = fopen(aliasfile_name, "r" );

    if( NULL == aliasfile )  //if the the protein alias file (.pal) file doesn't exist then
    {

        *num_volumes = 1; //there is only one volume

        char *indexfile_name = (char*)calloc(strlen(Database_copy) + strlen(".pin") + 1, sizeof(char));
        strcat(indexfile_name, Database_copy);
        strcat(indexfile_name, ".pin");
        FILE* indexfile = fopen(indexfile_name, "r" );//read index file

        if( NULL == indexfile)
        {
            fprintf(stderr, "Cannot open index file %s\n", indexfile_name);
            exit(1);
        }
        else
        {
            size_t result = 0;
            unsigned char* buffer = (char*)malloc(10000*sizeof(char));
            result = fread(buffer, sizeof(char), 10000, indexfile);//read the first 10000 bytes of the index file

            volume_lengths_temp = (Int4*) realloc( *volume_lengths, 1*sizeof(Int4) );
            if( NULL != volume_lengths_temp )
                *volume_lengths = volume_lengths_temp;
            else
            {
                printf("Error (re)allocating during index file reading\n");
                exit(1);
            }

            (*volume_lengths)[0] = ExtractVolumeSequences(buffer);

            free(buffer);
        }

        free(indexfile_name);
        fclose(indexfile);

    }
    else     //else read the lines line by line to find the "DBLIST" line
    {

        char* line_read = (char *) malloc( sizeof(char) );
        size_t nbytes = 1;
        Int4 bytes_read = getline(&line_read, &nbytes, aliasfile);//Read firs line

        char * pch;
        *num_volumes = 0;
        int k = 0;
        for(k = 0; k < 1000; ++k) //read the first 1000 lines of the alias file
        {

            if( feof(aliasfile) ) //at end of file break
                break;

            pch = strtok (line_read," "); //extract the first word of the line
            int diff = strncmp(pch, "DBLIST", 6); //compare the word with "DBLIST"
            if(diff)  //if it is not "DBLIST" read the next line and continue
            {

                bytes_read = getline(&line_read, &nbytes, aliasfile);//Read firs line which is a comment
                continue;

            }
            else     //if it is "DBLIST" read the rest of the volume names
            {

                pch = strtok (NULL, " "); //extract the first volume name
                while (pch != NULL)
                {

                    //check if the last character of pch is '\n' which means this is the last index file
                    if( pch[strlen(pch) - 1 ] == '\n')
                        pch[strlen(pch) - 1 ] = '\0'; //replace with the null character

                    char *indexfile_name = (char*)calloc(strlen(pch) + strlen(".pin") + 1, sizeof(char));
                    strcat(indexfile_name, pch);
                    strcat(indexfile_name, ".pin");
                    FILE* indexfile = fopen(indexfile_name, "r" );//read index file

                    if( NULL == indexfile)
                    {

                        fprintf(stderr, "Cannot open the index file %s\n", indexfile_name);
                        exit(1);

                    }
                    else
                    {
                        size_t result = 0;
                        unsigned char* buffer = (char*)malloc(10000*sizeof(char));
                        result = fread(buffer, sizeof(char), 1000, indexfile);//read the first 1000 bytes of the index file


                        volume_lengths_temp = (Int4*) realloc( *volume_lengths, ((*num_volumes)+1)*sizeof(Int4) );
                        if( NULL != volume_lengths_temp )
                            *volume_lengths = volume_lengths_temp;
                        else
                        {
                            printf("Error (re)allocating during index file reading\n");
                            exit(1);
                        }

                        (*volume_lengths)[*num_volumes] = ExtractVolumeSequences(buffer);

                        printf("%s: %d\n",indexfile_name, (*volume_lengths)[*num_volumes]);
                        free(buffer);

                    }

                    (*num_volumes)++;

                    free(indexfile_name);
                    fclose(indexfile);
                    pch = strtok (NULL, " "); //extract the next volume name

                }


                break; //after done with reading all volume names break
            }

        }

        free(line_read);
    }

    if( NULL != aliasfile )
        fclose(aliasfile);

    free(aliasfile_name);
    free(Database_copy);

}


/** Creates the database used by the H_BLAST
 @param Sequence_length_vector vector that holds the length of each sequence group [in]
 @param Group_number           the number of groups that the sequence database is split [in]
 @param Group_size             size of each group of sequences [in]
 @param stride                 that is used in createing the GPU database [in]
 @param seq_arc                data structure that holds information about the sequence database [in]
 @param num_sequence_to        last sequence of the database that is stored in the GPU database [in]
*/
void CreateH_BLAST_database(const BlastSeqSrc* seq_src,
                               const BlastGPUOptions* gpu_options)
{


    Int4 i = 0, j = 0;
    Int4 stride = 4;
    Int4 Group_size = (gpu_options->num_blocksx) * (gpu_options->num_threadsx);
    char percentage = PERCENTAGE;

    BlastSeqSrcIterator* itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    BlastSeqSrcGetSeqArg seq_arg;
    memset((void*) &seq_arg, 0, sizeof(seq_arg));

    /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
       sequences are retrieved in ncbistdaa/ncbi2na encodings respectively. */
    seq_arg.encoding = eBlastEncodingProtein;

    char* Database_name = BlastSeqSrcGetName(seq_src);
    char* GPU_Database_name = NULL;
    char* GPU_Database_name_temp;

    char* H_BLAST_database = NULL;
    char* H_BLAST_database_temp;
    Int4* Sequence_length_vector = NULL;
    Int4* Sequence_length_vector1 = NULL;
    Int4* Sequence_length_vector_temp;
    char* digits = (char*) calloc(10, sizeof(char));//10 digits are enough to cover up to 999999999 volumes


    char* GPU_information_name = (char*) calloc(strlen(Database_name) + strlen(".hDbI") + 1,sizeof(char));
    strcat(GPU_information_name, Database_name);
    strcat(GPU_information_name, ".hDbI");

    FILE* GPU_information_file = fopen( GPU_information_name, "wb");
    if ( NULL == GPU_information_file )
        printf("GPU information file cannot be opened for writing\n");

    Int4* volume_lengths = NULL;
    Int4 num_volumes = 0;
    ReadH_BLAST_volinfo(seq_src, &num_volumes, &volume_lengths);
    printf("num_volumes = %d\n",num_volumes);
    Write_GPU_information_file_preliminary(GPU_information_file, gpu_options->num_blocksx,
                                           gpu_options->num_threadsx, percentage, num_volumes);

    //Calculate the total size of the rest of the CPU database
    unsigned long CPU_rest_bytes = 0;

    Int4 volume = 0, break_limit = 0;
    for(volume = 0; volume < num_volumes; ++volume)
    {
        CPU_rest_bytes = 0;
        Int4 num_sequences_to = (99.9 * volume_lengths[volume]) / 100;
        num_sequences_to = (num_sequences_to / 32)*32;
        Int4 store_limit = break_limit + num_sequences_to;
        Int4 Group_number = (num_sequences_to / Group_size ) + (( (num_sequences_to % Group_size ) == 0)?0:1);

        Sequence_length_vector_temp = (Int4*) realloc( Sequence_length_vector, Group_number*sizeof(Int4) );
        if( NULL != Sequence_length_vector_temp )
            Sequence_length_vector = Sequence_length_vector_temp;
        else
        {
            printf("Error (re)allocating memory for the GPU database creation\n");
            exit(1);
        }

        Int4 vector_counter = 0;
        unsigned long H_BLAST_database_bytes = 0;
        Int4 oid = break_limit + Group_size - 1; //initialized oid to the last sequence of the first group
        if( oid >= store_limit )
            oid = store_limit - 1;

        //Initializes the Sequence_length_vector with the length of the longest sequence in each group
        //The length is converted to a multiple of WIDTH_MULTIPLE (defined in gpu_cpu_common.h)
        int previous_sequence_length = 0;
        while( oid <= store_limit - 1 )
        {

            Int4 sequence_length = BlastSeqSrcGetSeqLen( seq_src, (void*) &oid );
            if( previous_sequence_length > sequence_length )
            {
                fprintf(stderr,"ERROR: the input database is not sorted\n");
                fprintf(stderr,"Sort the database by using the option \"-sort_volumes\" with \"makeblastdb.\"\n");
                fprintf(stderr,"Exiting...\n");
                exit(1);
            }
            previous_sequence_length = sequence_length;

            Round2Multiple( &sequence_length );
            Sequence_length_vector[ vector_counter ] = sequence_length;
            if( GPU_VERBOSE )
                printf("Group %d starts at sequence %d and has length = %d\n",
                       vector_counter, oid, Sequence_length_vector[vector_counter]);

            ++vector_counter;

            oid += Group_size;

        }


        //If the last iteration of the previous do-while didn't cover all sequences (i.e. the Group_size is not a multiple of num_sequences)
        //the remaining sequences form the last group, and the last element of Sequence_length_vector is initialized
        if( oid != store_limit - 1 )
        {
            oid = store_limit - 1;
            Int4 sequence_length = BlastSeqSrcGetSeqLen( seq_src, (void*) &oid );
            /* check if the volume is sorted */
            if( previous_sequence_length > sequence_length )
            {
                fprintf(stderr,"ERROR: the input database is not sorted\n");
                fprintf(stderr,"Sort the database by using the option \"-sort_volumes\" with \"makeblastdb.\"\n");
                fprintf(stderr,"Exiting...\n");
                exit(1);
            }
            Round2Multiple( &sequence_length );
            Sequence_length_vector[ vector_counter ] = sequence_length;
            if( GPU_VERBOSE )
                printf("Group %d starts at sequence %d and has length = %d\n",
                       vector_counter, oid, Sequence_length_vector[vector_counter]);
        }


        for(i = (num_sequences_to / Group_size ) -1 ; i >=0 ; i--)
            H_BLAST_database_bytes += Sequence_length_vector[i]*Group_size*sizeof(char);
        H_BLAST_database_bytes += ( (num_sequences_to % Group_size ) >0 )?
                                     (( (num_sequences_to % Group_size )/(gpu_options->num_threadsx) + ((num_sequences_to%(gpu_options->num_threadsx) >0)?1:0) )*gpu_options->num_threadsx*Sequence_length_vector[Group_number-1]):0;




        int total_offset = 0;
        for(i = 0; i < Group_number; ++i)
            total_offset += Sequence_length_vector[i]+40;
        total_offset -=40;
        if(total_offset>= (32767) )
        {
            printf("Error: over 32k!\n");
        }
        else
        {
            printf("total_offse is %d\n", total_offset);
        }

        int total_ok_group = 0, total_ok_sub_group=0;



        unsigned int aux_group_size = gpu_options->num_threadsx * 16;
        printf("volumne %d\n", volume);
        printf("The rest sequence count = %d\n", volume_lengths[volume]-num_sequences_to+1);
        printf("The rest sequence group = %d\n", (volume_lengths[volume]-num_sequences_to+1)/aux_group_size);


        int off_count = 0, imb_group_count=0;
        unsigned long min_length=0, max_length=0;
        double imb_ratio = 0;

        unsigned long imb_mask = (1<<10)-1;


        Int4 count = 0;
        oid = store_limit;
        for(count = num_sequences_to; count<volume_lengths[volume]; count++)
        {
            CPU_rest_bytes += BlastSeqSrcGetSeqLen( seq_src, (void*) &oid ) * sizeof(char);
            if(off_count==0)
                min_length = BlastSeqSrcGetSeqLen( seq_src, (void*) &oid );

            off_count++;
            if((off_count & imb_mask) == 0)
            {

                if(off_count>0)
                {
                    max_length = BlastSeqSrcGetSeqLen( seq_src, (void*) &oid );
                    if(min_length==0)
                    {
                        printf("Error, min_length is zero!\n");
                        exit(-1);
                    }
                    imb_ratio = (double)max_length / min_length;
                    if(imb_ratio>1.1)
                        imb_group_count++;
                    min_length = max_length;
                    if((max_length+total_offset+40)<32767)
                        total_ok_sub_group ++;
                }
            }

            if(((off_count/16384)*16384) == off_count)
            {
                total_offset += max_length+40;
                if(total_offset< 32767)
                    total_ok_group ++;
                else
                    printf("x");
            }


            oid++;

        }

        int actual_total_length = 0;


        if( H_BLAST_database_bytes > UINT4_MAX )
        {
            printf("Memory requirements for GPU database larger than 4 GB. Exiting...\n");
            exit(0);
        }

        H_BLAST_database_temp = (char*) realloc( H_BLAST_database, H_BLAST_database_bytes*sizeof(char) );
        if( NULL != H_BLAST_database_temp )
            H_BLAST_database = H_BLAST_database_temp;
        else
        {
            printf("Error (re)allocating memory for the GPU database creation\n");
            exit(1);
        }

        memset(H_BLAST_database, 0, H_BLAST_database_bytes*sizeof(char));

        int length = 0;
        int group_member = 0;
        int counter = 0, sequence_id_in_group = 0, AminoAcids = 0, base = 0;
        int sequence_id_in_block=0, block_id=0, block_base=0;
        break_limit += volume_lengths[volume];
        /* iterate over all subject sequences of this volume */
        previous_sequence_length = 0;
        while ( ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) )
        {
            if (seq_arg.oid == BLAST_SEQSRC_ERROR)
                break;
            if (BlastSeqSrcGetSequence(seq_src,  &seq_arg) < 0)
                continue;


            if( seq_arg.oid < store_limit )
            {

                if( (counter % Group_size) == 0 )
                {
                    base += length*Group_size;
                    length = Sequence_length_vector[ group_member ];
                    group_member++;
                    sequence_id_in_group = 0;
                    block_id = 0;
                    block_base = 0;
                }

                if( (counter % (gpu_options->num_threadsx)) == 0 )
                {
                    sequence_id_in_block = 0;
                    if( (counter % Group_size) != 0 )
                    {
                        block_id++;
                        block_base += gpu_options->num_threadsx * length;
                    }
                }

                AminoAcids = 0;
                int sequence_length = BlastSeqSrcGetSeqLen( seq_src, (void*) &seq_arg.oid );
                if( previous_sequence_length > sequence_length )
                {
                    fprintf(stderr,"ERROR: the input database is not sorted\n");
                    fprintf(stderr,"Sort the database by using the option \"-sort_volumes\" with \"makeblastdb.\"\n");
                    fprintf(stderr,"Exiting...\n");
                    exit(1);
                }
                previous_sequence_length = sequence_length;

                for(i = 0; i < sequence_length/stride + 1; ++i)
                {
                    for(j = 0; j < stride; ++j)
                    {
                        (H_BLAST_database)[ base + block_base + (sequence_id_in_block + i*gpu_options->num_threadsx)*stride + j] = (seq_arg.seq->sequence)[AminoAcids];
                        AminoAcids++;
                        if( AminoAcids >= sequence_length )
                            goto done;
                    }
                }


done:
                actual_total_length += sequence_length;

                counter++;
                sequence_id_in_group++;
                sequence_id_in_block++;

            }

            BlastSeqSrcReleaseSequence(seq_src, &seq_arg);

            if( seq_arg.oid == (break_limit-1) )
                break;

        } /* while loop */


        //Create the name of the gpu database
        if( 1 == num_volumes )
        {
            int GPU_Database_name_bytes = strlen(Database_name) + strlen(".hDb") + 1;
            GPU_Database_name_temp = (char*) realloc(GPU_Database_name, GPU_Database_name_bytes*sizeof(char));
            if( NULL != GPU_Database_name_temp )
            {
                GPU_Database_name = GPU_Database_name_temp;
                memset(GPU_Database_name, 0, GPU_Database_name_bytes);
            }
            else
            {
                printf("Error (re)allocating memory for the GPU database creation\n");
                exit(1);
            }
            strcat(GPU_Database_name, Database_name);
        }
        else     //if there are more than one volumes then extra characters are needed for the volume number
        {
            sprintf(digits,".%02d",volume);
            Int4 GPU_Database_name_bytes = strlen(Database_name) + strlen(digits) + strlen(".hDb") + 1;
            GPU_Database_name_temp = (char*) realloc(GPU_Database_name, GPU_Database_name_bytes*sizeof(char));
            if( NULL != GPU_Database_name_temp )
            {
                GPU_Database_name = GPU_Database_name_temp;
                memset(GPU_Database_name, 0, GPU_Database_name_bytes);
            }
            else
            {
                printf("Error (re)allocating memory for the GPU database creation\n");
                exit(1);
            }
            strcat(GPU_Database_name, Database_name);
            strcat(GPU_Database_name, digits);
        }

        strcat(GPU_Database_name, ".hDb");

        FILE* Database_file = fopen( GPU_Database_name,"wb");
        if ( NULL == Database_file )
            printf("GPU Database file cannot be opened for writing\n");

        if( GPU_VERBOSE )
            printf("GPU BLASTP database size = %ld bytes\n", H_BLAST_database_bytes);

        fwrite(H_BLAST_database, sizeof(char), H_BLAST_database_bytes, Database_file);
        fclose(Database_file);

        Write_GPU_information_file(GPU_information_file,
                                   GPU_Database_name,
                                   H_BLAST_database_bytes,
                                   volume_lengths[volume],
                                   num_sequences_to,
                                   store_limit,
                                   break_limit,
                                   Group_number,
                                   Sequence_length_vector,
                                   CPU_rest_bytes);


        printf("Done with creating the GPU Database file (%s)\n", GPU_Database_name);

    }

    fclose(GPU_information_file);
    free(GPU_Database_name);
    free(H_BLAST_database);
    free(Sequence_length_vector);
    free(digits);

}

void H_BLAST_Execute(
    const BlastSeqSrc* seq_src,
    const LookupTableWrap* lookup_wrap,
    const BlastCoreAuxStruct* aux_struct,
    const BLAST_SequenceBlk* query,
    const BlastQueryInfo* query_info,
    const BlastInitialWordParameters* word_params,
    const BlastGPUOptions* gpu_options,
    const Int4 h_Hits_bytes,
    const Int4 h_GPUBlastInitHitList_bytes,
    const char* GPU_volume_name,
    const unsigned long h_H_BLAST_database_bytes,
    Int4* h_Sequence_length_vector,
    const Int4 Group_number,
    double *GPU_elapsed,
    double *FillDatabase_elapsed,
    int volume_id)
{

    struct timespec GPU_start, GPU_end, CPU_start, CPU_end;

    clock_gettime(CLOCK_MONOTONIC, &CPU_start);

    if(GPU_VERBOSE)
        fprintf(stderr,"threadId = 0: Reading the GPU database...\n");

    Read_H_BLAST_GPU_DB(
        GPU_volume_name,
        h_H_BLAST_database_bytes,
        volume_id);
    if(GPU_VERBOSE)
        fprintf(stderr,"threadId = 0: Done reading the GPU database\n");

    clock_gettime(CLOCK_MONOTONIC, &CPU_end);
    *FillDatabase_elapsed = (CPU_end.tv_sec - CPU_start.tv_sec) + (double)(CPU_end.tv_nsec - CPU_start.tv_nsec)/ BILLION;

    clock_gettime(CLOCK_MONOTONIC, &GPU_start);


    H_BLAST_pre(
        query, query_info, lookup_wrap, aux_struct->ewp, word_params,
        gpu_options, h_Sequence_length_vector,
        h_H_BLAST_database_bytes,
        Group_number,
        h_Hits_bytes, h_GPUBlastInitHitList_bytes
    );

    clock_gettime(CLOCK_MONOTONIC, &GPU_end);
    *GPU_elapsed = (double) (GPU_end.tv_sec - GPU_start.tv_sec) + (double)(GPU_end.tv_nsec - GPU_start.tv_nsec) / BILLION;

}


Int4
GPU_BLAST_PreliminarySearchEngine(EBlastProgramType program_number,
                                  BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
                                  const BlastSeqSrc* seq_src, BlastGapAlignStruct* gap_align,
                                  BlastScoringParameters* score_params,
                                  LookupTableWrap* lookup_wrap,
                                  const BlastInitialWordOptions* word_options,
                                  BlastExtensionParameters* ext_params,
                                  BlastHitSavingParameters* hit_params,
                                  BlastEffectiveLengthsParameters* eff_len_params,
                                  const PSIBlastOptions* psi_options,
                                  const BlastDatabaseOptions* db_options,
                                  const BlastGPUOptions* gpu_options,
                                  BlastHSPStream* hsp_stream, BlastDiagnostics* diagnostics,
                                  TInterruptFnPtr interrupt_search, SBlastProgress* progress_info)
{
    BlastCoreAuxStruct* aux_struct = NULL;
    BlastHSPList* hsp_list = NULL;
    BlastSeqSrcGetSeqArg seq_arg;
    BlastSeqSrcGetSeqArg gpu_seq_arg;
    Int2 status = 0;
    Int8 db_length = 0;
    const BlastScoringOptions* score_options = score_params->options;
    const BlastHitSavingOptions* hit_options = hit_params->options;
    const BlastExtensionOptions* ext_options = ext_params->options;
    BlastInitialWordParameters* word_params = NULL;
    Boolean gapped_calculation = score_options->gapped_calculation;
    BlastScoreBlk* sbp = gap_align->sbp;
    BlastSeqSrcIterator* itr;
    BlastSeqSrcIterator* gpu_itr;


    memset((void*) &seq_arg, 0, sizeof(seq_arg));
    memset((void*) &gpu_seq_arg, 0, sizeof(gpu_seq_arg));


    const Boolean kNucleotide = (program_number == eBlastTypeBlastn ||
                                 program_number == eBlastTypePhiBlastn);

    BlastInitialWordParametersNew(program_number, word_options,
                                  hit_params, lookup_wrap, sbp, query_info,
                                  BlastSeqSrcGetAvgSeqLen(seq_src), &word_params);

    if ((status =
                s_BlastSetUpAuxStructures(seq_src, lookup_wrap, word_params,
                                          ext_options, hit_options, query, &aux_struct)) != 0)
        return status;

    pthread_mutex_lock(&thread_counter_mutex);
    Int4 threadId = thread_counter++;

    if( 0 == threadId )
    {
        int init_return = pthread_barrier_init(&barr, NULL, gpu_options->cpu_threads);

        Read_H_BLAST_runtime_options();
        init_max_and_current_gpu_walltime();
        lut_sorting((BlastAaLookupTable*) lookup_wrap->lut);


        int num_contexts_invoked = query_info->last_context - query_info->first_context +1;
        int avg_query_length = 0;
        int q = 0;

        avg_query_length  = query->length / num_contexts_invoked;

        length_per_query=avg_query_length;


    }
    pthread_mutex_unlock(&thread_counter_mutex);


    itr = BlastSeqSrcIteratorNewEx(subject_chunck_size);
    gpu_itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/4500,1));

    if ( (-1 == gpu_options->method) )
    {
        if( 0 == threadId )  //only one threads creates the database
        {
            CreateH_BLAST_database(seq_src, gpu_options);
        }
        exit(0);
    }
    else
    {
        /// Set no. of gpu cards for alignment
        if( 0 == threadId )
        {
            set_first_gpu_id(gpu_options->use_gpu -1);

            if( (0 == gpu_options->method) && (gpu_options->use_gpu>1))
            {
                printf("Warning: the \'-gpu 0\' is not compatible with the option \'-method 0\'!\n");
                printf("The option \'-method 0\' has changed to \'-method 1\'!\n");
                SetNumOfGpuCards(1);
            }
            else SetNumOfGpuCards( gpu_options->method);
        }

    }

    pthread_barrier_wait(&barr);

    struct timespec CPU_start, CPU_end, Total_start, Total_end;

    double GPU_elapsed = 0, CPU_elapsed = 0, Total_elapsed = 0, FillDatabase_elapsed = 0;

    Int4 num_blocksx, num_threadsx, num_volumes = 0, percentage = 0;
    FILE* GPU_information_file = NULL;

    init_max_and_current_cpu_walltime();

    float local_gpu_db_used_ratio = Gpu_db_used_ratio_in_byte;

    Read_GPU_DB_information_file_preliminary(&GPU_information_file, seq_src,
                                          &num_blocksx, &num_threadsx,
                                          &percentage, &num_volumes);

    BlastGPUOptions* gpu_options_local = (BlastGPUOptions *) calloc(2, sizeof(BlastGPUOptions));
    gpu_options_local[0].use_gpu = gpu_options->use_gpu;
    gpu_options_local[1].use_gpu = gpu_options->use_gpu;
    gpu_options_local[0].num_blocksx = num_blocksx;
    gpu_options_local[1].num_blocksx = num_blocksx;
    gpu_options_local[0].num_threadsx = num_threadsx;
    gpu_options_local[1].num_threadsx = num_threadsx;
    gpu_options_local[0].num_sequences_to = percentage;
    gpu_options_local[1].num_sequences_to = percentage;
    gpu_options_local[0].method = gpu_options->method;
    gpu_options_local[1].method = gpu_options->method;
    char* GPU_volume_name = NULL;
    Int4* Sequence_length_vector = NULL;

    //Initialize global variables
    if( 0 == threadId )
    {
        h_Hits = NULL;
        h_GPUBlastInitHitList = NULL;
        sufficient_memory[0] = FALSE;
        sufficient_memory[1] = FALSE;
    }

    Int4 volume = 0;

    float local_gpu_db_used_ratio_array[2] = {0.0, 0.0};
    local_gpu_db_used_ratio_array[0] = Gpu_db_used_ratio_in_byte;
    long int local_successful_extension_array[2] = {0,0};

    unsigned long h_H_BLAST_database_bytes[2];
    Int4 gpu_limit[2]= {0,0}, gpu_limit_temp;
    Int4 num_sequences_to[2]= {0,0};
    Int4 Group_number[2]= {0,0};
    Int4 gpu_NextChunkOID[2]= {0,0};
    Int4 volume_length[2]= {0,0};
    Int4 break_limit[2]= {0,0};

    /// for @gpu_db_used_ratio adjustment
    unsigned long local_db_volume_size_in_byte[3];

    init_db_volume_size_in_byte(local_db_volume_size_in_byte);

    Int4 first_oid[4] = {0,0,0,0};

    int SubjectIsTranslated_flag = Blast_SubjectIsTranslated(program_number);

    ///  Set MmapMemory
    if( (0 == threadId) && (gpu_options->cpu_threads > 1) )
    {
        BlastSeqSrcSetMmapMemory(seq_src, FALSE);

    }
    else if ( (0 == threadId)  )
    {
        BlastSeqSrcSetMmapMemory(seq_src, TRUE);
    }

    pthread_barrier_wait(&barr);


    /// Setup the init. gpu kernel selection
    if (0 == threadId)
    {
        if(query->length <= 1000)
            set_gpu_kernel_option(0);
        else
            set_gpu_kernel_option(1);
    }


    for(volume = -1; volume < num_volumes; ++volume)
    {
        Int4 h_Hits_bytes;
        Int4 h_GPUBlastInitHitList_bytes;
        Int4 break_limit_temp = 0;
        double H_BLAST_get_data_elapsed, freeing_memory, waiting_time1;

        clock_gettime(CLOCK_MONOTONIC, &Total_start);

        if(volume >-1 && sufficient_memory[CURRENT(volume)] )
        {
            clock_gettime(CLOCK_MONOTONIC, &CPU_start);

            pthread_barrier_wait(&barr);

            if ( (volume>0) && (0==threadId ) )
            {
                int q=0;
                double max_cpu_wall_time_in_previous_round = 0;

                for(q=0; q<gpu_options->cpu_threads; q++)
                {
                    if(current_cpu_walltime[q]>max_cpu_wall_time_in_previous_round)
                        max_cpu_wall_time_in_previous_round = current_cpu_walltime[q];
                }

                max_cpu_walltime[(volume-1)%3] = max_cpu_wall_time_in_previous_round;
            }

            clock_gettime(CLOCK_MONOTONIC, &CPU_end );
            waiting_time1 = ( CPU_end.tv_sec - CPU_start.tv_sec ) + (double)( CPU_end.tv_nsec - CPU_start.tv_nsec ) / BILLION;


            if( 0 == threadId )
            {
                pthread_mutex_lock(&thread_counter_mutex);

                get_GPU_result( &h_Hits,
                                h_Hits_bytes,
                                &h_GPUBlastInitHitList,
                                h_GPUBlastInitHitList_bytes,
                                num_sequences_to[CURRENT(volume)],
                                &H_BLAST_get_data_elapsed);

                clock_gettime(CLOCK_MONOTONIC, &CPU_start);
                H_BLAST_free_memory( );
                clock_gettime(CLOCK_MONOTONIC, &CPU_end );
                freeing_memory = ( CPU_end.tv_sec - CPU_start.tv_sec ) + (double)( CPU_end.tv_nsec - CPU_start.tv_nsec ) / BILLION;

                pthread_mutex_unlock(&thread_counter_mutex);
            }
        }
        pthread_barrier_wait(&barr);

        if(volume>=0)
            local_successful_extension_array[CURRENT(volume)]=get_successful_extension_in_volume();

        int gpu_db_group_size = gpu_options_local[0].num_blocksx * gpu_options_local[0].num_threadsx;
        if(volume<(num_volumes-1))
        {
            h_H_BLAST_database_bytes[NEXT(volume)] = 0;
            num_sequences_to[NEXT(volume)] = 0;
            gpu_limit[NEXT(volume)] = 0, gpu_limit_temp = 0;
            Group_number[NEXT(volume)] = 0;
            gpu_NextChunkOID[NEXT(volume)] = gpu_NextChunkOID[CURRENT(volume)] + volume_length[CURRENT(volume)];

            Read_GPU_DB_information_file(GPU_information_file,
                                      &GPU_volume_name,
                                      &(h_H_BLAST_database_bytes[NEXT(volume)]), &(volume_length[NEXT(volume)]),
                                      &(num_sequences_to[NEXT(volume)]), &gpu_limit_temp,
                                      &break_limit_temp, &(Group_number[NEXT(volume)]),
                                      &(local_db_volume_size_in_byte[(volume+1)%3]),
                                      &Sequence_length_vector);

            if ( (volume>0) && (flag_fixed_gpu_db_used_ratio!=1) )
            {
                local_gpu_db_used_ratio=get_new_gpudb_ratio(volume, local_gpu_db_used_ratio,  gpu_options->cpu_threads,
                                        get_num_of_gpu_cards_used(), local_db_volume_size_in_byte,
                                        local_gpu_db_used_ratio_array, local_successful_extension_array,
                                        gapped_calculation, query->length);
            }

            get_GPU_DB_information(
                &(h_H_BLAST_database_bytes[NEXT(volume)]), &(volume_length[NEXT(volume)]),
                &(num_sequences_to[NEXT(volume)]), &gpu_limit_temp,
                &break_limit_temp, &(Group_number[NEXT(volume)]),
                &Sequence_length_vector,
                gpu_options->cpu_threads,
                gpu_options_local[0].num_threadsx,
                gpu_db_group_size,
                &local_gpu_db_used_ratio);

            local_gpu_db_used_ratio_array[NEXT(volume)] = local_gpu_db_used_ratio;

            //if( (0 == threadId) && (volume>0) ) printf("vol = %d, mnew ratio = %lg\n", volume, local_gpu_db_used_ratio);

            gpu_options_local[NEXT(volume)].num_sequences_to = num_sequences_to[NEXT(volume)];
            gpu_limit[NEXT(volume)] = break_limit[CURRENT(volume)] + num_sequences_to[NEXT(volume)];
            break_limit[NEXT(volume)] = break_limit_temp;

            if(volume<(num_volumes-2))
                first_oid[(volume+2)&0x3] = break_limit[NEXT(volume)];

            h_Hits_bytes = num_sequences_to[NEXT(volume)] * H_HITS_SIZE_BYTE  * sizeof(Int4);

            h_GPUBlastInitHitList_bytes =  num_sequences_to[NEXT(volume)] * get_MAX_SUCCESSFUL_HITS_PER_SEQUENCE() * sizeof(GPUBlastInitHitList);
            if ( get_num_of_gpu_cards_used() > 1 ) h_GPUBlastInitHitList_bytes *= 1.5;

            sufficient_memory[NEXT(volume)] = H_BLAST_check_memory(lookup_wrap, aux_struct->ewp,
                                              &(gpu_options_local[NEXT(volume)]),
                                              h_H_BLAST_database_bytes[NEXT(volume)],
                                              h_Hits_bytes, h_GPUBlastInitHitList_bytes,
                                              Group_number[NEXT(volume)], query_info->num_queries,
                                              query->length, mem_ecc_flag);

            if( (0 == threadId) && !sufficient_memory[NEXT(volume)] )
            {
                fprintf(stderr,"WARNING: Not enough GPU global memory to process volume No. %02d of the database. Continuing without the GPU...\n",volume);
                fprintf(stderr,"         Consider splitting the input database in volumes with smaller size by using the option \"-max_file_size <String>\" (e.g. -max_file_sz 500MB).\n");
                fprintf(stderr,"         Consider using fewer GPU blocks and then fewer GPU thread when formating the database with (e.g. \"-gpu_blocks 256\" and/or \"-gpu_threads 192\") to reduce the GPU global memory requiremenets.\n");
            }


            if( sufficient_memory[NEXT(volume)] && (0 == threadId) )
            {
                if(BRIEF_INFO)
                    printf("Checking DB volume %d\n", volume+1);

                H_BLAST_Execute(seq_src, lookup_wrap, aux_struct, query, query_info, word_params,
                                   &(gpu_options_local[NEXT(volume)]),  h_Hits_bytes, h_GPUBlastInitHitList_bytes,
                                   GPU_volume_name, h_H_BLAST_database_bytes[NEXT(volume)],
                                   Sequence_length_vector, Group_number[NEXT(volume)],
                                   &GPU_elapsed, &FillDatabase_elapsed, volume+1);
            }
        }


        if(volume>-1)
        {

            pthread_barrier_wait(&barr);

            /// RPS processing
            /* remember the current search state */
            if (progress_info)
                progress_info->stage = ePrelimSearch;

            /* For RPS BLAST, there is no loop over subject sequences, so the preliminary
               search engine is done in a separate function. */
            if (Blast_ProgramIsRpsBlast(program_number))
            {
                status =
                    s_RPSPreliminarySearchEngine(program_number, query, query_info,
                                                 seq_src, score_params, lookup_wrap, aux_struct, word_params,
                                                 ext_params, gap_align, hit_params, hsp_stream, diagnostics,
                                                 interrupt_search, progress_info,
                                                 &(gpu_options_local[CURRENT(volume)]), h_Hits, h_GPUBlastInitHitList );
                word_params = BlastInitialWordParametersFree(word_params);
                s_BlastCoreAuxStructFree(aux_struct);
                return status;
            }

            /// CPU post-processing preperation
            /* Update the parameters for linking HSPs, if necessary. */
            BlastLinkHSPParametersUpdate(word_params, hit_params, gapped_calculation);

            /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
               sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
            seq_arg.encoding = eBlastEncodingProtein;

            db_length = BlastSeqSrcGetTotLen(seq_src);


            //if there is not enough GPU memory use only the CPU and start alignment from the first database sequence
            num_sequences_to[CURRENT(volume)] = sufficient_memory[CURRENT(volume)] ? num_sequences_to[CURRENT(volume)] : 0;
            gpu_limit[CURRENT(volume)] = sufficient_memory[CURRENT(volume)] ? gpu_limit[CURRENT(volume)] : (break_limit[CURRENT(volume)] - volume_length[CURRENT(volume)]);

            Int4 print_point = 0;

            clock_gettime(CLOCK_MONOTONIC, &CPU_start);
            while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
                    != BLAST_SEQSRC_EOF)
            {

                if (seq_arg.oid == BLAST_SEQSRC_ERROR)
                    break;

                if( sufficient_memory[CURRENT(volume)]  &&  ( seq_arg.oid < gpu_limit[CURRENT(volume)] ) && (h_Hits[(seq_arg.oid - first_oid[volume&0x3])*4+1] == 0) )
                    continue;

                if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0)
                    continue;

                if (db_length == 0)
                {
                    /* This is not a database search, hence need to recalculate and save
                       the effective search spaces and length adjustments for all
                       queries based on the length of the current single subject
                       sequence. */
                    if ((status = BLAST_OneSubjectUpdateParameters(program_number,
                                  seq_arg.seq->length, score_options, query_info,
                                  sbp, hit_params, word_params,
                                  eff_len_params)) != 0)
                        return status;
                }


                /* Calculate cutoff scores for linking HSPs. Do this only for
                   ungapped protein searches and ungapped translated
                   searches. */
                if (hit_params->link_hsp_params && !kNucleotide && !gapped_calculation)
                {
                    CalculateLinkHSPCutoffs(program_number, query_info, sbp,
                                            hit_params->link_hsp_params, word_params, db_length,
                                            seq_arg.seq->length);
                }



                /// For CPU DB processing
                if (seq_arg.oid >= gpu_limit[CURRENT(volume)] && ( SubjectIsTranslated_flag >0 ) )
                {
                    Uint4 i;
                    for (i = 0; i < seq_arg.seq->num_seq_ranges; i++)
                    {
                        seq_arg.seq->seq_ranges[i].left =
                            CONV_NUCL2PROT_COORDINATES(seq_arg.seq->seq_ranges[i].left);
                        seq_arg.seq->seq_ranges[i].right =
                            CONV_NUCL2PROT_COORDINATES(seq_arg.seq->seq_ranges[i].right);
                    }


                    if (seq_arg.seq->gen_code_string == NULL)
                    {
                        seq_arg.seq->gen_code_string =
                            GenCodeSingletonFind(db_options->genetic_code);
                    }
                    ASSERT(seq_arg.seq->gen_code_string);
                }

                /// CPU post-processing
                if( sufficient_memory[CURRENT(volume)]  &&  ( seq_arg.oid < gpu_limit[CURRENT(volume)] ) )
                {
                    gpu_options_local[CURRENT(volume)].use_gpu = TRUE;

                    status =
                        s_BlastSearchEngineCore(program_number, query, query_info,
                                                seq_arg.seq, lookup_wrap, gap_align, score_params, word_params,
                                                ext_params, hit_params, db_options, diagnostics, aux_struct,
                                                &hsp_list, interrupt_search, progress_info,
                                                &(gpu_options_local[CURRENT(volume)]), h_Hits, h_GPUBlastInitHitList, gpu_NextChunkOID[CURRENT(volume)]);


                }
                /// CPU DB processing
                else if(seq_arg.oid >= gpu_limit[CURRENT(volume)])
                {
                    gpu_options_local[CURRENT(volume)].use_gpu = FALSE;

                    status =
                        s_BlastSearchEngineCore(program_number, query, query_info,
                                                seq_arg.seq, lookup_wrap, gap_align, score_params, word_params,
                                                ext_params, hit_params, db_options, diagnostics, aux_struct,
                                                &hsp_list, interrupt_search, progress_info,
                                                &(gpu_options_local[CURRENT(volume)]), h_Hits, h_GPUBlastInitHitList, 0);
                }


                if (status)
                {
                    break;
                }

                if (hsp_list && hsp_list->hspcnt > 0)
                {
                    if (!gapped_calculation)
                    {
                        /* The following must be performed for any ungapped
                           search with a nucleotide database. */
                        status =
                            Blast_HSPListReevaluateUngapped(
                                program_number, hsp_list, query,
                                seq_arg.seq, word_params, hit_params,
                                query_info, sbp, score_params, seq_src,
                                seq_arg.seq->gen_code_string);

                        if (status)
                        {
                            BlastSeqSrcReleaseSequence(seq_src, &seq_arg);
                            return status;
                        }
                        /* Relink HSPs if sum statistics is used, because scores might
                         * have changed after reevaluation with ambiguities, and there
                         * will be no traceback stage where relinking is done normally.
                         * If sum statistics are not used, just recalculate e-values.
                         */
                        if (hit_params->link_hsp_params)
                        {
                            status =
                                BLAST_LinkHsps(program_number, hsp_list, query_info,
                                               seq_arg.seq->length, sbp,
                                               hit_params->link_hsp_params,
                                               gapped_calculation);
                        }
                        else
                        {
                            Blast_HSPListGetEvalues(program_number,query_info,
                                                    seq_arg.seq->length, hsp_list,
                                                    gapped_calculation,
                                                    0, sbp, 0, 1.0);
                        }
                        /* Use score threshold rather than evalue if
                            * matrix_only_scoring is used.  -RMH-
                            */
                        if ( sbp->matrix_only_scoring )
                        {
                            status = Blast_HSPListReapByRawScore(hsp_list,
                                                                 hit_params->options);
                        }
                        else
                        {
                            status = Blast_HSPListReapByEvalue(hsp_list,
                                                               hit_params->options);
                        }

                        /* Calculate and fill the bit scores, since there will be no
                           traceback stage where this can be done. */
                        Blast_HSPListGetBitScores(hsp_list, gapped_calculation, sbp);
                    }

                    /* Save the results. */
                    status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
                    if (status != 0)
                        break;

                    if (hit_params->low_score)
                    {
                        int query_index;
                        for (query_index=0; query_index<hsp_stream->results->num_queries; query_index++)
                            if (hsp_stream->results->hitlist_array[query_index] && hsp_stream->results->hitlist_array[query_index]->heapified)
                                hit_params->low_score[query_index] =
                                    MAX(hit_params->low_score[query_index],
                                        hit_params->options->low_score_perc*(hsp_stream->results->hitlist_array[query_index]->low_score));
                    }

                } //if (hsp_list && hsp_list->hspcnt > 0)

                BlastSeqSrcReleaseSequence(seq_src, &seq_arg);

                /* check for interrupt */
                if (interrupt_search && (*interrupt_search)(progress_info) == TRUE)
                {
                    status = BLASTERR_INTERRUPTED;
                    break;
                }

                if( seq_arg.oid >= (break_limit[CURRENT(volume)]-1) )
                    break;
            }

            clock_gettime(CLOCK_MONOTONIC, &CPU_end);
            CPU_elapsed = (CPU_end.tv_sec - CPU_start.tv_sec) + (double)(CPU_end.tv_nsec - CPU_start.tv_nsec) / BILLION;


            clock_gettime(CLOCK_MONOTONIC, &CPU_start);

            //wait until both threads are done with BLASTP_GetGPUHits() before deallocating h_Hits and h_GPUBlastInitHitList since
            //if one threadId=0 deallocates these arrays while other threads are still working on them a segmentation fault will occur
            if( sufficient_memory[CURRENT(volume)] )
                pthread_barrier_wait(&barr);

            clock_gettime(CLOCK_MONOTONIC, &CPU_end );
            double waiting_time3 = ( CPU_end.tv_sec - CPU_start.tv_sec ) + (double)( CPU_end.tv_nsec - CPU_start.tv_nsec ) / BILLION;

            clock_gettime(CLOCK_MONOTONIC, &Total_end);
            Total_elapsed = (Total_end.tv_sec - Total_start.tv_sec) + (double)(Total_end.tv_nsec - Total_start.tv_nsec)/ BILLION;
            current_cpu_walltime[threadId] = CPU_elapsed;

        }

    }

    if( sufficient_memory[NEXT(volume)] && (0 == threadId) )
    {
        free(h_Hits);
        free(h_GPUBlastInitHitList);
    }

    pthread_mutex_lock(&thread_counter_mutex);
    thread_counter--;
    pthread_mutex_unlock(&thread_counter_mutex);

    free(Sequence_length_vector);
    free(GPU_volume_name);
    free(gpu_options_local);
    fclose(GPU_information_file);

    if(0==threadId)
        all_gpu_quit();


    hsp_list = Blast_HSPListFree(hsp_list);  /* in case we were interrupted */
    BlastSequenceBlkFree(seq_arg.seq);
    itr = BlastSeqSrcIteratorFree(itr);
    gpu_itr = BlastSeqSrcIteratorFree(gpu_itr);

    /* Fill the cutoff values in the diagnostics structure */
    if (diagnostics && diagnostics->cutoffs)
    {
        s_FillReturnCutoffsInfo(diagnostics->cutoffs, score_params, word_params,
                                ext_params, hit_params);
    }

    word_params = BlastInitialWordParametersFree(word_params);
    s_BlastCoreAuxStructFree(aux_struct);

    return status;
}

