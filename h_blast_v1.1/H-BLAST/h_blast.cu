#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "h_blast.h"
#include "h_blast_kernel.cu"
#include <algo/blast/core/gpu_cpu_common.h>
#include "h_blast_kernel.h"

#include <algo/blast/core/blast_aalookup.h>

#include <algorithm>
#include <vector>


// includes for mmap
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

// include for mmap mamagement
#include <algo/blast/core/mmap_storage_manager.h>


using namespace std;

#define CUDA_DEVICE_MAXNUMBER 100
/// The macro "NO_CELL_LEVELS" shows the statistic
/// function will report at most 2^(NO_CELL_LEVELS-1)
/// results can be stored within a subject.
#define NO_CELL_LEVELS 7


typedef struct
{
    int GPU_device_id;
    char* h_db_access_entry;
    int* h_db_Sequence_Maxlength_vector;
    Int4 number_of_Groups;
    Int4 num_seq_to;
    unsigned long h_H_BLAST_database_bytes;
    Int4 h_Hits_offset;
    Int4 h_GPUBlastInitHitList_offset;
    Int4 h_Hits_bytes;
    Int4 h_GPUBlastInitHitList_bytes;
    Int4* h_Hits;
    GPUBlastInitHitList* h_GPUBlastInitHitList;
    Int4* d_Hits;
    GPUBlastInitHitList* d_GPUBlastInitHitList;
    Int4* d_Database2Dpadded;
    Int4* d_RepeatedSubstitutionMatrix;
    Int4* d_RepeatedSequence_length_vector;
    PV_ARRAY_TYPE* d_RepeatedPV;
    AaLookupSmallboneCell* d_ThickBackbone;
    Uint2* d_overflow;
    short* d_RepeatedDiag_array;

    cudaStream_t stream_id;
    cudaEvent_t event_id, start;

    int* d_hash_table;
    int* d_one_hit_buffer;

} GPU_Stream_Block;

static unsigned int MAX_SUCCESSFUL_HITS_PER_SEQUENCE=19;

/// The option shows which gpu kernel is used.
/// '0' for the first one, while '1' for the other.
static int gpu_kernel_option = 0;


/// The option shows the first gpu device id used in this
/// process. We assume that all gpus are in exculsive mode.
static unsigned int first_gpu_id = 0;

static int prefixsum_hash_table4[16] = {0,585, 4680, 1169, 8768, 1161, 5256, 1745, 12800, 1097, 5192, 1681, 9280, 1673, 5768, 2257};

int get_processor_no()
{
    int ax =0;
    FILE *fp;
    char *line_read = (char *) malloc( sizeof(char) );
    size_t nbytes = 1;
    int bytes_read;
    fp=popen("cat /proc/cpuinfo | grep \"physical id\" | sort | uniq | wc -l", "r");
    bytes_read = getline(&line_read, &nbytes, fp);

    if(bytes_read>0) ax = atoi(line_read);
    else ax=1;

    free(line_read);
    pclose(fp);
    return ax;
}


int get_core_count_per_processor()
{
    int ax =0;
    FILE *fp;
    char *line_read = (char *) malloc( sizeof(char) );
    size_t nbytes = 1;
    int bytes_read;
    fp=popen("grep \"cpu cores\" /proc/cpuinfo |sort -u |cut -d\":\" -f2", "r");

    bytes_read = getline(&line_read, &nbytes, fp);

    if(bytes_read>0) ax = atoi(line_read);
    else ax=1;

    free(line_read);
    pclose(fp);
    return ax;
}

void set_first_gpu_id(unsigned int id)
{
    first_gpu_id = id;
}


void set_gpu_kernel_option(int flag)
{
    if( (flag==0) || (flag==1))
        gpu_kernel_option = flag;
    return;
}

int static get_gpu_kernel_option()
{
    return gpu_kernel_option;
}

static long int successful_extension_in_volume = 0;

long int get_successful_extension_in_volume ()
{
    return  successful_extension_in_volume;
}

void static small_bobbo_sorting(unsigned short* short_array, int array_count )
{
    int i,j;
    Uint2 tmp=0;
    for(i=0; i<(array_count-1); i++)
        for(j=0; j<(array_count-i); j++)
        {
            if(short_array[j]<short_array[j-1])
            {
                tmp = short_array[j];
                short_array[j] = short_array[j-1];
                short_array[j-1] = tmp;
            }
        }
    return ;
}


void lut_sorting(BlastAaLookupTable *lookup)
{


    AaLookupSmallboneCell* h_ThickBackbone = (AaLookupSmallboneCell*) lookup->thick_backbone;

    unsigned short* h_overflow = (unsigned short*)lookup->overflow;

    unsigned long int h_overflow_count = lookup->overflow_size;

    vector<unsigned short> overflow_vector (h_overflow, h_overflow+h_overflow_count);

    int i;

    for(i=0; i<lookup->backbone_size; i++)
    {
        if(h_ThickBackbone[i].num_used > AA_HITS_PER_CELL)
            sort(overflow_vector.begin()+ h_ThickBackbone[i].payload.overflow_cursor,  overflow_vector.begin() + h_ThickBackbone[i].payload.overflow_cursor + h_ThickBackbone[i].num_used  );
        else if(h_ThickBackbone[i].num_used > 0 )
            small_bobbo_sorting( h_ThickBackbone[i].payload.entries,  h_ThickBackbone[i].num_used);
    }

    for(i=0; i<h_overflow_count ; i++)
    {
        h_overflow[i] = overflow_vector[i];
    }

    return;

}


static GPU_Stream_Block stream_blocks[CUDA_DEVICE_MAXNUMBER];
static double current_gpu_walltime[CUDA_DEVICE_MAXNUMBER];
static double max_gpu_walltime=0.0;

double get_gpu_walltime();

void set_max_gpu_walltime();

void init_max_and_current_gpu_walltime();

int get_num_of_gpu_cards_used();

static int num_of_gpu_cards_used;

Int4* local_h_Hits=NULL;
GPUBlastInitHitList* local_h_GPUBlastInitHitList=NULL;

char* h_Database2Dpadded=NULL;
unsigned long h_H_BLAST_database_bytes_one_volume = 0;
int file_id_for_gpu_db=0;

void gpu_set_device(const int device_id);
void gpu_get_device_count(int * count);


double get_gpu_walltime()
{
    return max_gpu_walltime;
}

int get_num_of_gpu_cards_used()
{
    return num_of_gpu_cards_used;
}
void set_max_gpu_walltime()
{
    int i;
    double max_record=0;
    for(i=0; i<num_of_gpu_cards_used; i++)
    {
        if(max_record<current_gpu_walltime[i])
            max_record = current_gpu_walltime[i];
    }
    max_gpu_walltime= max_record;
}

void init_max_and_current_gpu_walltime()
{
    int i;
    for(i=0; i<num_of_gpu_cards_used; i++)
        current_gpu_walltime[i]=0;

    max_gpu_walltime=0;
}

void set_MAX_SUCCESSFUL_HITS_PER_SEQUENCE(const unsigned int val)
{
    MAX_SUCCESSFUL_HITS_PER_SEQUENCE = val;
}

unsigned int get_MAX_SUCCESSFUL_HITS_PER_SEQUENCE()
{
    return MAX_SUCCESSFUL_HITS_PER_SEQUENCE;
}

void SetNumOfGpuCards(const int expected_num_of_cards)
{
    gpu_get_device_count(&num_of_gpu_cards_used);
    if(num_of_gpu_cards_used>CUDA_DEVICE_MAXNUMBER)
    {
        printf("Warning: The no. of avalible gpu cards are more then %d.\nNow, uses %d cards only!\nTo use all cards, please modify the constant CUDA_DEVICE_MAXNUMBER\n",
               CUDA_DEVICE_MAXNUMBER, num_of_gpu_cards_used);
        num_of_gpu_cards_used = CUDA_DEVICE_MAXNUMBER;
    }

    if(0==expected_num_of_cards)
    {
        if(first_gpu_id>1)
        {
            printf("Warning: There is no gpu card avalible, uses cpus only!\n" );
            num_of_gpu_cards_used = 0;
        }

    }
    else
    {
        int delta = num_of_gpu_cards_used -first_gpu_id;
        if(delta>0)
        {
            if(delta>=expected_num_of_cards)
            {
                num_of_gpu_cards_used = expected_num_of_cards;
            }
            else
            {
                printf("Warning: There is not enough gpu cards avalible, uses %d cards only!\n", delta );
                num_of_gpu_cards_used = delta;
            }
        }
        else
        {
            printf("Warning: There is no gpu card avalible!\n" );
            num_of_gpu_cards_used = 0;
        }
    }
}


void gpu_quit()
{
    cudaDeviceReset();
}


void gpu_abnormally_exit(int exit_code)
{
    gpu_quit();
    exit(exit_code);
}

void all_gpu_quit()
{
    int gpu_device_count=num_of_gpu_cards_used, i;

    for(i=0; i<gpu_device_count; i++)
    {
        gpu_set_device(stream_blocks[i].GPU_device_id);
        cudaDeviceReset();
    }
}

void all_gpu_abnormally_exit(int exit_code)
{
    all_gpu_quit();
    exit(exit_code);
}


void gpu_set_cache_config(enum cudaFuncCache cacheConfig)
{
    cudaError_t return_code = cudaFuncSetCacheConfig(H_BLAST_kernelTwoHit2, cacheConfig);
    if(cudaSuccess  != return_code)
    {
        fprintf(stderr, "cudaFuncSetCacheConfig error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    return;
}

void gpu_memcpy_async(void * dst,
                      const void *  src,
                      size_t  count,
                      enum cudaMemcpyKind  kind,
                      cudaStream_t stream)
{
    cudaError_t return_code = cudaMemcpyAsync (dst, src, count, kind, stream);
    if(cudaSuccess  != return_code)
    {
        fprintf(stderr, "cudaMemcpyAsync error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    return;
}


void gpu_memcpy(void * dst,
                const void *  src,
                size_t  count,
                enum cudaMemcpyKind  kind
               )
{
    cudaError_t return_code = cudaMemcpy(dst, src, count, kind);
    if(cudaSuccess  != return_code)
    {
        fprintf(stderr, "cudaMemcpy error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    return;
}



void gpu_set_device(const int device_id)
{
    cudaError_t return_code = cudaSetDevice(device_id);
    if(cudaSuccess  != return_code)
    {

        fprintf(stderr, "cudaSetDevice error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    return;

}

void gpu_device_mem_free(void ** d_ptr)
{
    if(NULL==*d_ptr)
    {
        fprintf(stderr, "gpu_device_mem_free error: input pointer is NULL\n");
        all_gpu_abnormally_exit(-1);
    }
    cudaError_t return_code = cudaFree(*d_ptr);
    if(cudaSuccess != return_code)
    {

        fprintf(stderr, "cudaFree error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    *d_ptr = NULL;
    return;
}


void gpu_host_calloc(void ** ptr,
                     size_t  size,
                     unsigned int  flags=cudaHostAllocDefault
                    )
{
    cudaError_t return_code = cudaHostAlloc(ptr, size, flags);
    if(cudaSuccess != return_code)
    {
        fprintf(stderr, "cudaHostAlloc error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    memset(*ptr, 0, size);
    return;

}


void gpu_device_calloc(void ** d_ptr, size_t size)
{
    cudaError_t return_code = cudaMalloc(d_ptr, size);
    if(cudaSuccess != return_code)
    {
        fprintf(stderr, "cudaMalloc error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    cudaMemset(*d_ptr, 0, size);
    return;

}


void gpu_host_mem_free(void** h_ptr)
{

    if(NULL==*h_ptr)
    {
        fprintf(stderr, "gpu_host_mem_free error: input pointer is NULL\n");
        all_gpu_abnormally_exit(-1);
    }
    cudaError_t return_code = cudaFreeHost(*h_ptr);
    if( cudaSuccess !=return_code)
    {
        fprintf(stderr, "cudaFreeHost error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    *h_ptr = NULL;
    return ;
}

#ifndef CUDA55
void gpu_memcpy_to_symbol(const char *  symbol,
#else
void gpu_memcpy_to_symbol(const void *  symbol,
#endif
                          const void *  src,
                          size_t  count,
                          size_t  offset=0,
                          enum cudaMemcpyKind kind=cudaMemcpyHostToDevice)
{
    cudaError_t return_code = cudaMemcpyToSymbol(symbol, src, count, offset,kind);
    if(cudaSuccess  != return_code)
    {

        fprintf(stderr, "cudaMemcpyToSymbol error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    return;
}


void gpu_get_device_count(int * count)
{
    cudaError_t return_code = cudaGetDeviceCount(count);
    if(cudaSuccess  != return_code)
    {

        fprintf(stderr, "cudaGetDeviceCount error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    return;
}


void init_stream_obj(int cuda_device_no)
{
    stream_blocks[cuda_device_no].GPU_device_id = cuda_device_no + first_gpu_id;

    gpu_set_device(stream_blocks[cuda_device_no].GPU_device_id);

    cudaError_t return_code=cudaStreamCreate(&stream_blocks[cuda_device_no].stream_id);

    if( cudaSuccess !=return_code)
    {
        fprintf(stderr, "cudaStreamCreate error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }

    return_code=cudaEventCreate(&stream_blocks[cuda_device_no].event_id);

    if( cudaSuccess !=return_code)
    {
        fprintf(stderr, "cudaEventCreate error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }

    return_code=cudaEventCreate(&stream_blocks[cuda_device_no].start);

    if( cudaSuccess !=return_code)
    {
        fprintf(stderr, "cudaEventCreate error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }

    return;
}

void GPU_Host_Alloc_WC(void** h_ptr, const size_t h_size)
{
    cudaError_t return_code=cudaHostAlloc(h_ptr,  h_size, cudaHostAllocWriteCombined);

    if( cudaSuccess !=return_code)
    {
        fprintf(stderr, "cudaHostAlloc error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }
    return;
}


void gpu_stream_synchronize(int cuda_device_id)
{

    cudaError_t return_code;

    gpu_set_device(stream_blocks[cuda_device_id].GPU_device_id);

    return_code=cudaEventSynchronize(stream_blocks[cuda_device_id].event_id);
    if( cudaSuccess !=return_code)
    {
        fprintf(stderr, "cudaEventSynchronize error: %s\n", cudaGetErrorString(return_code));
        all_gpu_abnormally_exit(-1);
    }

    float kernel_time;
    cudaEventElapsedTime(&kernel_time, stream_blocks[cuda_device_id].start, stream_blocks[cuda_device_id].event_id);

#ifdef DEBUG
    printf("Kernel elapses %lg sec\n", kernel_time/1000);
#endif

    current_gpu_walltime[cuda_device_id] = kernel_time/1000;

    return;
}


void Clean_H_BLAST_database()
{
    /// return state
    int ret = -1;
    if(h_Database2Dpadded!=NULL)
    {
        ret = munmap(h_Database2Dpadded, h_H_BLAST_database_bytes_one_volume);
        if(-1 == ret)
        {
            fputs("munmap previous gpu db volume failed: \n", stderr);
            all_gpu_abnormally_exit (1);
        }
        close(file_id_for_gpu_db);
        h_H_BLAST_database_bytes_one_volume = 0;
        h_Database2Dpadded = NULL;
    }

}


void Read_H_BLAST_GPU_DB(
    const char* GPU_volume_name,
    const unsigned long h_H_BLAST_database_bytes,
    const int volume_id)
{
    if(is_mmap_allocated_by_id(volume_id)>0)
    {
        mmap_storage_manager tmp =  get_mmap_storage_manager_by_id(volume_id);

        h_Database2Dpadded = tmp.start_for_mem_file;
        h_H_BLAST_database_bytes_one_volume = h_H_BLAST_database_bytes;
        file_id_for_gpu_db = tmp.file_id_for_mem_file;

        if(tmp.size_of_mem_file_in_bytes < h_H_BLAST_database_bytes)
        {
            fputs("The length of gpu db is less than required\n",stderr);
            all_gpu_abnormally_exit (1);
        }

        return;
    }

    /// file handle
    int file_id=0;
    /// the pointer for using file content in memory
    char* start_for_mem_file=NULL;
    /// return state
    int ret = -1;
    /// file status, including the length of the file
    struct stat buf = {0};

    file_id=open(GPU_volume_name,O_RDONLY);

    ret = fstat(file_id, &buf);
    if(-1 == ret)
    {
        fputs("fstat failed:\n",stderr);
        //fputs("The file %s is not existed!\n",GPU_volume_name,stderr);
        all_gpu_abnormally_exit (1);
    }

    if(buf.st_size < h_H_BLAST_database_bytes)
    {
        fputs("The length of gpu db is less than required\n",stderr);
        all_gpu_abnormally_exit (1);
    }

    start_for_mem_file=(char*)mmap(NULL,buf.st_size,PROT_READ,MAP_SHARED|MAP_POPULATE,file_id,0);

    if(start_for_mem_file==MAP_FAILED)
    {
        fputs("File error: GPU database file could not be oppened\n",stderr);
        all_gpu_abnormally_exit (1);
    }

    h_Database2Dpadded = start_for_mem_file;
    h_H_BLAST_database_bytes_one_volume = h_H_BLAST_database_bytes;
    file_id_for_gpu_db = file_id;


    /// Add the mmap to the existed pool
    if(push_a_mmap_storage_manager_to_array_with_info(start_for_mem_file, buf.st_size, file_id) < 0)
    {
        fputs("Error: Can not add the mmap to the pool!\n",stderr);
        all_gpu_abnormally_exit (1);
    }
}


void setup_data_partitions(Int4 Group_number,
                           int * Sequence_length_vector,
                           char * h_H_BLAST_database,
                           int cuda_device_count,
                           const Int4 Group_size,
                           const Int4 num_seq_to,
                           const Int4 total_sequence_bytes,
                           const Int4 total_init_hsp_list_length)
{

    Int4 average_sequence_bytes=0, sequence_bytes_of_a_part=0,  group_ptr=0, last_group_ptr=0;
    Int4 last_total_database_bytes=0;
    Int4 i=0;

    if(cuda_device_count==0)
    {
        printf("Error: the cuda device count is 0!\n");
        exit(-1);
    }

    average_sequence_bytes=total_sequence_bytes*1.0/cuda_device_count;
    sequence_bytes_of_a_part = 0;

    last_group_ptr=-1;
    group_ptr=-1;
    last_total_database_bytes=0;

    if (1==cuda_device_count)
    {
        stream_blocks[0].h_db_Sequence_Maxlength_vector=Sequence_length_vector;
        stream_blocks[0].number_of_Groups=Group_number;
        stream_blocks[0].h_H_BLAST_database_bytes = total_sequence_bytes ;
        stream_blocks[0].h_Hits_offset =0;
        stream_blocks[0].h_GPUBlastInitHitList_offset = 0;
        stream_blocks[0].num_seq_to = num_seq_to ;
        stream_blocks[0].h_db_access_entry= h_H_BLAST_database;
    }
    else
    {

        char flag=0;

        for(i=0; i< cuda_device_count-1; i++)
        {
            stream_blocks[i].h_db_access_entry= &(h_H_BLAST_database[last_total_database_bytes]);

            while(group_ptr +1< Group_number-1)
            {
                sequence_bytes_of_a_part +=  Sequence_length_vector[group_ptr+1]*Group_size;

                if(sequence_bytes_of_a_part>average_sequence_bytes)
                {
                    stream_blocks[i].h_db_Sequence_Maxlength_vector=&(Sequence_length_vector[last_group_ptr+1]);
                    stream_blocks[i].number_of_Groups=group_ptr-last_group_ptr;
                    stream_blocks[i].num_seq_to = stream_blocks[i].number_of_Groups*Group_size;
                    if(group_ptr==last_group_ptr)
                    {
                        printf("Error: The size of a group of subject sequences is greater than the amount of average sequence bytes of a part! \
                               Please decrease the gpu thread size or the gpu block size!\n");
                        exit(-1);
                    }
                    else
                    {
                        stream_blocks[i].h_H_BLAST_database_bytes = sequence_bytes_of_a_part - Sequence_length_vector[group_ptr+1]*Group_size;
                    }

                    stream_blocks[i].h_Hits_offset = (last_group_ptr+1)*Group_size*H_HITS_SIZE_BYTE;
                    stream_blocks[i].h_GPUBlastInitHitList_offset = i*total_init_hsp_list_length/cuda_device_count;

                    last_total_database_bytes += stream_blocks[i].h_H_BLAST_database_bytes;
                    sequence_bytes_of_a_part = 0;

                    last_group_ptr=group_ptr;
                    break;
                }
                group_ptr ++;
            }

            if(group_ptr+1 == Group_number-1)
            {
                if(0==flag)
                {
                    stream_blocks[i].h_H_BLAST_database_bytes = total_sequence_bytes - last_total_database_bytes;

                    if((group_ptr==last_group_ptr) && (stream_blocks[i].h_H_BLAST_database_bytes > average_sequence_bytes))
                    {
                        printf("Error: The size of the last group of subject sequences is greater than the amount of average sequence bytes of a part! \
                               Please decrease the gpu thread size or the gpu block size!\n");
                        exit(-1);
                    }

                    else if(stream_blocks[i].h_H_BLAST_database_bytes > average_sequence_bytes)
                    {
                        flag = 1;
                        stream_blocks[i].h_H_BLAST_database_bytes = sequence_bytes_of_a_part;
                        stream_blocks[i].h_db_Sequence_Maxlength_vector=&(Sequence_length_vector[last_group_ptr+1]);
                        stream_blocks[i].number_of_Groups=group_ptr-last_group_ptr;
                        stream_blocks[i].h_Hits_offset = (last_group_ptr+1)*Group_size*H_HITS_SIZE_BYTE;
                        stream_blocks[i].h_GPUBlastInitHitList_offset = i*total_init_hsp_list_length/cuda_device_count;
                        stream_blocks[i].num_seq_to = stream_blocks[i].number_of_Groups*Group_size;
                        last_group_ptr=group_ptr;
                        last_total_database_bytes += stream_blocks[i].h_H_BLAST_database_bytes;
                        sequence_bytes_of_a_part=0;
                    }

                    else
                    {

                        stream_blocks[i].h_db_Sequence_Maxlength_vector=&(Sequence_length_vector[Group_number-1]);
                        stream_blocks[i].number_of_Groups=1;

                        stream_blocks[i].h_Hits_offset = (last_group_ptr+1)*Group_size*H_HITS_SIZE_BYTE;
                        stream_blocks[i].h_GPUBlastInitHitList_offset =  i*total_init_hsp_list_length/cuda_device_count;
                        stream_blocks[i].num_seq_to = num_seq_to - (last_group_ptr+1)*Group_size;

                        flag = 2;
                        group_ptr++;
                        last_total_database_bytes = total_sequence_bytes;
                        last_group_ptr = Group_number-1;
                    }
                }
                else if(1==flag)
                {

                    stream_blocks[i].h_db_Sequence_Maxlength_vector=&(Sequence_length_vector[last_group_ptr+1]);
                    stream_blocks[i].number_of_Groups=1;
                    stream_blocks[i].h_H_BLAST_database_bytes = total_sequence_bytes - last_total_database_bytes;
                    if(stream_blocks[i].h_H_BLAST_database_bytes > average_sequence_bytes)
                    {
                        printf("Error: The size of the last group of subject sequences is greater than the \n amount of average sequence bytes of a part! \n \
                               Please decrease the gpu thread size or the gpu block size!\n");
                        exit(-1);
                    }
                    stream_blocks[i].h_Hits_offset = (Group_number-1)*Group_size*H_HITS_SIZE_BYTE;
                    stream_blocks[i].h_GPUBlastInitHitList_offset = i*total_init_hsp_list_length/cuda_device_count;

                    stream_blocks[i].num_seq_to = num_seq_to - (Group_number-1)*Group_size;

                    flag=2;
                    group_ptr++;
                    last_total_database_bytes = total_sequence_bytes;
                    last_group_ptr = Group_number-1;
                }

            }

            if(group_ptr+1 == Group_number)
            {
                stream_blocks[i].h_db_Sequence_Maxlength_vector=Sequence_length_vector;
                stream_blocks[i].number_of_Groups=0;
                stream_blocks[i].h_H_BLAST_database_bytes = 0;
                stream_blocks[i].h_Hits_offset = 0;
                stream_blocks[i].h_GPUBlastInitHitList_offset = 0;
                stream_blocks[i].num_seq_to = 0;
            }

        }

        if(group_ptr+1 == Group_number)
        {
            stream_blocks[cuda_device_count-1].h_db_access_entry= h_H_BLAST_database;
            stream_blocks[cuda_device_count-1].h_db_Sequence_Maxlength_vector=Sequence_length_vector;
            stream_blocks[cuda_device_count-1].number_of_Groups=0;
            stream_blocks[cuda_device_count-1].h_H_BLAST_database_bytes = 0;
            stream_blocks[cuda_device_count-1].h_Hits_offset = 0;
            stream_blocks[cuda_device_count-1].h_GPUBlastInitHitList_offset = 0;
            stream_blocks[cuda_device_count-1].num_seq_to = 0;
        }
        else
        {
            stream_blocks[cuda_device_count-1].h_db_access_entry= &(h_H_BLAST_database[last_total_database_bytes]);
            stream_blocks[cuda_device_count-1].h_db_Sequence_Maxlength_vector=&(Sequence_length_vector[last_group_ptr+1]);
            stream_blocks[cuda_device_count-1].number_of_Groups=Group_number-last_group_ptr-1;
            stream_blocks[cuda_device_count-1].h_H_BLAST_database_bytes = total_sequence_bytes - last_total_database_bytes;


            stream_blocks[cuda_device_count-1].h_Hits_offset = (last_group_ptr+1)*Group_size*H_HITS_SIZE_BYTE;
            stream_blocks[cuda_device_count-1].h_GPUBlastInitHitList_offset = (cuda_device_count-1)*total_init_hsp_list_length/cuda_device_count;
            stream_blocks[cuda_device_count-1].num_seq_to = num_seq_to - (last_group_ptr+1)*Group_size;
        }
    }

    for(i=0; i<cuda_device_count; i++)
    {
        stream_blocks[i].h_Hits_bytes =   stream_blocks[i].num_seq_to * H_HITS_SIZE_BYTE * sizeof(Int4);
        stream_blocks[i].h_GPUBlastInitHitList_bytes = total_init_hsp_list_length * sizeof(GPUBlastInitHitList) / cuda_device_count;
        stream_blocks[i].h_Hits = NULL;
        stream_blocks[i].h_GPUBlastInitHitList = NULL;
    }

    Int4 DB_bytes = 0, total_seq_num_to =0, total_db_group=0;

    for(i=0; i<cuda_device_count; i++)
    {
        DB_bytes += stream_blocks[i].h_H_BLAST_database_bytes;
        total_seq_num_to += stream_blocks[i].num_seq_to;
        total_db_group += stream_blocks[i].number_of_Groups;
    }

    if(DB_bytes != total_sequence_bytes)
        printf("total db bytes is wrong!\n");

    if(total_seq_num_to != num_seq_to)
        printf("total_seq_num_to is wrong!\n");

    if(total_db_group != Group_number)
        printf("total_db_group is wrong!\n");


    for(i=0; i<cuda_device_count; i++)
    {
        if(stream_blocks[i].num_seq_to==0) stream_blocks[i].GPU_device_id = -1;
    }
    return;
}


void  H_BLAST_pre(
    const BLAST_SequenceBlk* query,
    const BlastQueryInfo* query_info,
    const LookupTableWrap* lookup_wrap,
    const Blast_ExtendWord* ewp,
    const BlastInitialWordParameters* word_params,
    const BlastGPUOptions* gpu_options,
    Int4* Sequence_length_vector,
    const int h_Database2Dpadded_bytes,
    const Int4 Group_number,
    const Int4 h_Hits_bytes,
    const Int4 h_GPUBlastInitHitList_bytes
)
{

    if(local_h_Hits!=NULL) gpu_host_mem_free((void**)&local_h_Hits);
    if(local_h_GPUBlastInitHitList!=NULL) gpu_host_mem_free((void**)&local_h_GPUBlastInitHitList);

    gpu_host_calloc((void**)&local_h_Hits, h_Hits_bytes, cudaHostAllocPortable);
    gpu_host_calloc((void**)&local_h_GPUBlastInitHitList, h_GPUBlastInitHitList_bytes, cudaHostAllocPortable);

    int gpu_device_count=num_of_gpu_cards_used;


    setup_data_partitions(Group_number,
                          Sequence_length_vector,
                          h_Database2Dpadded,
                          gpu_device_count,
                          gpu_options->num_blocksx * gpu_options->num_threadsx,
                          gpu_options->num_sequences_to,
                          h_Database2Dpadded_bytes,
                          h_GPUBlastInitHitList_bytes/sizeof(GPUBlastInitHitList));

    for(int cuda_device_id=0; cuda_device_id<gpu_device_count; cuda_device_id++)
    {

        init_stream_obj(cuda_device_id);

        //gpu_set_cache_config(cudaFuncCachePreferShared);
        H_BLAST(cuda_device_id,
                   query, query_info,
                   lookup_wrap,
                   ewp, word_params,
                   gpu_options,
                   stream_blocks[cuda_device_id].h_db_Sequence_Maxlength_vector
                  );
    }
}

void  H_BLAST(const int cuda_device_id,
                 const BLAST_SequenceBlk* query, const BlastQueryInfo* query_info,
                 const LookupTableWrap* lookup_wrap,
                 const Blast_ExtendWord* ewp, const BlastInitialWordParameters* word_params,
                 const BlastGPUOptions* gpu_options, const Int4* Sequence_length_vector
                )
{

    const Int4 num_threadsy = 1, num_blocksy = 1;
    const Int4 num_blocksx = gpu_options->num_blocksx, num_threadsx = gpu_options->num_threadsx,  num_sequences = stream_blocks[cuda_device_id].num_seq_to;

    if( GPU_VERBOSE )
        fprintf(stderr,"threadId = 0: Number of processed sequences by the GPU = %5d\n", num_sequences);

    int Group_number = (num_sequences / (num_threadsx*num_blocksx) ) + (( (num_sequences % (num_threadsx*num_blocksx) ) == 0)?0:1);

    int Sequence_length_vector_stride = ( Group_number / PITCH + ((Group_number % PITCH)?1:0) ) * PITCH;

    BlastAaLookupTable *lookup = (BlastAaLookupTable*) lookup_wrap->lut;

    PV_ARRAY_TYPE *h_pv = lookup->pv;
    int pv_length = (lookup->backbone_size >> PV_ARRAY_BTS) + 1;

    Round2Multiple(&pv_length);

    unsigned long int h_RepeatedPV_bytes = num_blocksx * pv_length * sizeof(PV_ARRAY_TYPE);

    PV_ARRAY_TYPE* h_RepeatedPV = (PV_ARRAY_TYPE*) calloc( h_RepeatedPV_bytes, sizeof(char));

    MultipleCopyPV(h_RepeatedPV, h_pv, pv_length, num_blocksx);

    AaLookupSmallboneCell* h_ThickBackbone = (AaLookupSmallboneCell*) lookup->thick_backbone;

    Uint2* h_overflow = (Uint2*)lookup->overflow;

    unsigned long int h_ThickBackbone_bytes = (lookup->backbone_size+1) * sizeof(AaLookupSmallboneCell);
    unsigned long int h_overflow_bytes = (lookup->overflow_size+1) * sizeof(Uint2);
    unsigned int diag_array_length = ewp->diag_table->diag_array_length;

    unsigned long long h_RepeatedDiag_array_bytes = (unsigned long long)diag_array_length * num_blocksx*num_threadsx*sizeof(short);
    if(0!=get_gpu_kernel_option())
        h_RepeatedDiag_array_bytes /=4;

    int diag_mask = ewp->diag_table->diag_mask;
    int window = ewp->diag_table->window;

    int Query_length = query->length;
    unsigned char* h_Query = query->sequence;

    int num_contexts_invoked = query_info->last_context - query_info->first_context +1;

    short* h_x_dropoff = (short*)malloc(num_contexts_invoked * sizeof(short) );
    short* h_cutoff_score = (short*)malloc(num_contexts_invoked * sizeof(short) );
    short* h_query_offset = (short*)malloc((num_contexts_invoked+1) * sizeof(short) );

    for(int i = 0; i < num_contexts_invoked; ++i)
    {
        h_x_dropoff[i] = (short)word_params->cutoffs[i].x_dropoff;
        h_cutoff_score[i] = (short)word_params->cutoffs[i].cutoff_score;
        h_query_offset[i] = (short)query_info->contexts[i].query_offset;
    }
    h_query_offset[num_contexts_invoked] = (short)query->length;

    char* h_SubstitutionMatrix = (char*) calloc( SUBSTITUTION_MATRIX_LENGTH * sizeof(char), sizeof(char) );

    unsigned int h_RepeatedSubstitutionMatrix_bytes = SUBSTITUTION_MATRIX_LENGTH * num_blocksx * sizeof(char);
    char* h_RepeatedSubstitutionMatrix = (char*) calloc ( h_RepeatedSubstitutionMatrix_bytes, sizeof(char) );


    ReadSubstitutionMatrix( h_SubstitutionMatrix, h_RepeatedSubstitutionMatrix, SUBSTITUTION_MATRIX_LENGTH, num_blocksx );

    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_Database2Dpadded), stream_blocks[cuda_device_id].h_H_BLAST_database_bytes + 2* num_threadsx * sizeof(int) ) ;
    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_RepeatedSubstitutionMatrix), h_RepeatedSubstitutionMatrix_bytes ) ;
    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_Hits), stream_blocks[cuda_device_id].h_Hits_bytes ) ;
    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_RepeatedPV), h_RepeatedPV_bytes ) ;
    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_ThickBackbone), h_ThickBackbone_bytes ) ;
    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_overflow), h_overflow_bytes ) ;
    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_RepeatedDiag_array), h_RepeatedDiag_array_bytes ) ;
    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_GPUBlastInitHitList), stream_blocks[cuda_device_id].h_GPUBlastInitHitList_bytes ) ;
    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_one_hit_buffer), ONE_HIT_BUFFER_SIZE*num_blocksx * num_threadsx * sizeof(int)/4 ) ;

    gpu_memcpy( stream_blocks[cuda_device_id].d_RepeatedSubstitutionMatrix, h_RepeatedSubstitutionMatrix, h_RepeatedSubstitutionMatrix_bytes, cudaMemcpyHostToDevice) ;
    gpu_memcpy( stream_blocks[cuda_device_id].d_RepeatedPV, h_RepeatedPV, h_RepeatedPV_bytes, cudaMemcpyHostToDevice) ;
    if(h_ThickBackbone != NULL)
        gpu_memcpy( stream_blocks[cuda_device_id].d_ThickBackbone, h_ThickBackbone, h_ThickBackbone_bytes, cudaMemcpyHostToDevice) ;
    if(h_overflow != NULL)
        gpu_memcpy( stream_blocks[cuda_device_id].d_overflow, h_overflow, h_overflow_bytes, cudaMemcpyHostToDevice) ;

    dim3  grid( num_blocksx, num_blocksy, 1);
    dim3  threads( num_threadsx, num_threadsy, 1);
    int  total_shared_mem_size_requested = 0, total_const_mem_size_requested=0;

    total_shared_mem_size_requested = pv_length*sizeof(PV_ARRAY_TYPE) + SUBSTITUTION_MATRIX_LENGTH*sizeof(char) + 5 * num_threadsx * 4 * sizeof(char) + num_threadsx*sizeof(int)+16*sizeof(int);
    total_const_mem_size_requested = Query_length + (num_contexts_invoked+1)*sizeof(int)*3 + Sequence_length_vector_stride*sizeof(int);

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, cuda_device_id);
    if(total_shared_mem_size_requested > deviceProp.sharedMemPerBlock)
    {
        fprintf(stderr,"ERROR: Shared memory requirements (%d bytes per block) bigger than available memory (%d bytes)!!\n",
                (int)total_shared_mem_size_requested, (int)deviceProp.sharedMemPerBlock);
        all_gpu_abnormally_exit(-1);
    }


    if(total_const_mem_size_requested > deviceProp.totalConstMem)
    {
        fprintf(stderr,"ERROR: Constant memory requirements (%d bytes per block) bigger than available memory (%d bytes)!!\n",
                (int)total_const_mem_size_requested, (int)deviceProp.totalConstMem);
        all_gpu_abnormally_exit(-1);
    }

    if((num_threadsx/32)*32 != num_threadsx)
    {
        fprintf(stderr,"ERROR: The no. of threads must be a multipy of 32!\n");
        all_gpu_abnormally_exit(-1);
    }

    int total_subject_length_offset_per_gpu_thread = 0;
    for(int i=0; i<Group_number; i++)
    {
        total_subject_length_offset_per_gpu_thread += Sequence_length_vector[i]+window;
    }

    if(total_subject_length_offset_per_gpu_thread > (1<<15) )
    {
        fprintf(stderr, "ERROR: the diag. array can not handle this query group!\n");
        all_gpu_abnormally_exit(-1);
    }

    int hash_subject_length = 0;
    for(int i=0; i<Group_number; i++)
    {
        hash_subject_length += Sequence_length_vector[i]+ewp->diag_table->window;
    }

    if(hash_subject_length>(1<<15))
    {
        fprintf(stderr,"ERROR: the  diag. array is not big enough, please make the size of each DB volumn smaller!\n") ;
        all_gpu_abnormally_exit(-1);
    }

    if(num_contexts_invoked>899)
    {
        fprintf(stderr,"ERROR: the arraies for query_info are not big enough, please processing less queries!\n") ;
        all_gpu_abnormally_exit(-1);
    }

#ifndef CUDA55
    gpu_memcpy_to_symbol( d_Query_const, h_Query, Query_length) ;
    gpu_memcpy_to_symbol( (char*)d_x_dropoff_const, h_x_dropoff, num_contexts_invoked*sizeof(short));
    gpu_memcpy_to_symbol( (char*)d_cutoff_score_const, h_cutoff_score, num_contexts_invoked*sizeof(short));
    gpu_memcpy_to_symbol( (char*)d_query_offset_const, h_query_offset, (num_contexts_invoked+1)*sizeof(short));
    gpu_memcpy_to_symbol( (char*)d_length_vector_const, Sequence_length_vector, Sequence_length_vector_stride*sizeof(int));
    gpu_memcpy_to_symbol( (char*)(&d_Database2Dpadded_const), &(stream_blocks[cuda_device_id].d_Database2Dpadded), sizeof(int*));
    gpu_memcpy_to_symbol( (char*)(&d_diag_array_const), &(stream_blocks[cuda_device_id].d_RepeatedDiag_array), sizeof(short*));
    gpu_memcpy_to_symbol( (char*)(&d_GPUBlastInitHitList_const), &(stream_blocks[cuda_device_id].d_GPUBlastInitHitList), sizeof(GPUBlastInitHitList*));
    gpu_memcpy_to_symbol( (char*)(&hit_list_base_const), &(stream_blocks[cuda_device_id].h_GPUBlastInitHitList_offset), sizeof(int));
    gpu_memcpy_to_symbol( (char*)(&diag_offset_const), &(window), sizeof(int));
#else
    gpu_memcpy_to_symbol( d_Query_const, h_Query, Query_length) ;
    gpu_memcpy_to_symbol( d_x_dropoff_const, h_x_dropoff, num_contexts_invoked*sizeof(short));
    gpu_memcpy_to_symbol( d_cutoff_score_const, h_cutoff_score, num_contexts_invoked*sizeof(short));
    gpu_memcpy_to_symbol( d_query_offset_const, h_query_offset, (num_contexts_invoked+1)*sizeof(short));
    gpu_memcpy_to_symbol( d_length_vector_const, Sequence_length_vector, Sequence_length_vector_stride*sizeof(int));
    gpu_memcpy_to_symbol( (&d_Database2Dpadded_const), &(stream_blocks[cuda_device_id].d_Database2Dpadded), sizeof(int*));
    gpu_memcpy_to_symbol( (&d_diag_array_const), &(stream_blocks[cuda_device_id].d_RepeatedDiag_array), sizeof(short*));
    gpu_memcpy_to_symbol( (&d_GPUBlastInitHitList_const), &(stream_blocks[cuda_device_id].d_GPUBlastInitHitList), sizeof(GPUBlastInitHitList*));
    gpu_memcpy_to_symbol( (&hit_list_base_const), &(stream_blocks[cuda_device_id].h_GPUBlastInitHitList_offset), sizeof(int));
    gpu_memcpy_to_symbol( (&diag_offset_const), &(window), sizeof(int));
#endif

    gpu_memcpy( & (stream_blocks[cuda_device_id].d_Database2Dpadded[ num_threadsx]), stream_blocks[cuda_device_id].h_db_access_entry,
                stream_blocks[cuda_device_id].h_H_BLAST_database_bytes,
                cudaMemcpyHostToDevice) ;

    gpu_device_calloc( (void**) &(stream_blocks[cuda_device_id].d_hash_table),16*sizeof(int)) ;
    gpu_memcpy( stream_blocks[cuda_device_id].d_hash_table, prefixsum_hash_table4, 16*sizeof(int),  cudaMemcpyHostToDevice);

    cudaEventRecord(stream_blocks[cuda_device_id].start,stream_blocks[cuda_device_id].stream_id);

    if(NULL == stream_blocks[cuda_device_id].h_db_access_entry)
    {
        printf("Error, the h_db entry is NULL, and the card id is %d\n", cuda_device_id);
        gpu_abnormally_exit(-1);
    }

    if(0==get_gpu_kernel_option())
    {
        H_BLAST_kernelTwoHit1<<< grid, threads,  total_shared_mem_size_requested, stream_blocks[cuda_device_id].stream_id >>>(
            stream_blocks[cuda_device_id].d_RepeatedSubstitutionMatrix,
            SUBSTITUTION_MATRIX_LENGTH,
            stream_blocks[cuda_device_id].d_RepeatedPV,
            stream_blocks[cuda_device_id].d_ThickBackbone,
            stream_blocks[cuda_device_id].d_overflow,
            Sequence_length_vector_stride,
            num_sequences, pv_length,
            Query_length, diag_mask, diag_array_length,
            num_contexts_invoked,
            lookup->word_length,
            lookup->mask, window,
            stream_blocks[cuda_device_id].d_Hits,
            stream_blocks[cuda_device_id].h_GPUBlastInitHitList_bytes/sizeof(GPUBlastInitHitList) );
    }
    else
    {
        H_BLAST_kernelTwoHit2<<< grid, threads,  total_shared_mem_size_requested, stream_blocks[cuda_device_id].stream_id >>>(
            stream_blocks[cuda_device_id].d_RepeatedSubstitutionMatrix,
            SUBSTITUTION_MATRIX_LENGTH,
            stream_blocks[cuda_device_id].d_RepeatedPV,
            stream_blocks[cuda_device_id].d_ThickBackbone,
            stream_blocks[cuda_device_id].d_overflow,
            Sequence_length_vector_stride,
            num_sequences, pv_length,
            Query_length, diag_mask, diag_array_length,
            num_contexts_invoked,
            lookup->word_length,
            lookup->mask, window,
            stream_blocks[cuda_device_id].d_Hits,
            stream_blocks[cuda_device_id].h_GPUBlastInitHitList_bytes/sizeof(GPUBlastInitHitList),
            stream_blocks[cuda_device_id].d_hash_table,
            stream_blocks[cuda_device_id].d_one_hit_buffer
        );
    }

    cudaEventRecord(stream_blocks[cuda_device_id].event_id, stream_blocks[cuda_device_id].stream_id);

    free(h_x_dropoff);
    free(h_cutoff_score);
    free(h_query_offset);
    free(h_SubstitutionMatrix);

    free(h_RepeatedPV);
    free(h_RepeatedSubstitutionMatrix);

}


void memory_defrag(
    const Int4 h_Hits_bytes,
    const Int4 h_GPUBlastInitHitList_bytes,
    const Int4 num_sequences)
{
    int i=0, current_hit_id=0, total_sucessful_extention_for_subject=0;
    GPUBlastInitHitList* GPUBlastInitHitList_tmp=
        (GPUBlastInitHitList*)malloc(h_GPUBlastInitHitList_bytes*sizeof(char));

    GPUBlastInitHitList* current_hit_list_in_this_subject=NULL;

    if((NULL==local_h_Hits) || (NULL==local_h_GPUBlastInitHitList))
    {
        printf("Error: the memory address of local_h_Hits or local_h_GPUBlastInitHitList is unaccessable!\n");
        gpu_abnormally_exit(-1);
    }
    if(NULL==GPUBlastInitHitList_tmp)
    {
        printf("Error: there is no enough memory space to varialbe GPUBlastInitHitList_tmp!\n");
        gpu_abnormally_exit(-1);
    }

    for(i=0; i<num_sequences; i++)
    {
        if( (-1 != local_h_Hits[i*4+3]) && (local_h_Hits[i*4+1]>0)  )
        {
            total_sucessful_extention_for_subject = local_h_Hits[i*4+1];
            current_hit_list_in_this_subject = & (local_h_GPUBlastInitHitList[local_h_Hits[i*4+3]]);
            local_h_Hits[i*4+3] = current_hit_id;

            if( 0==(i&0x10) )
            {
                memcpy(  GPUBlastInitHitList_tmp+current_hit_id,current_hit_list_in_this_subject,total_sucessful_extention_for_subject*sizeof(GPUBlastInitHitList));
            }
            else
            {
                memcpy(  GPUBlastInitHitList_tmp+current_hit_id,  current_hit_list_in_this_subject-total_sucessful_extention_for_subject+1, total_sucessful_extention_for_subject*sizeof(GPUBlastInitHitList));
            }

            current_hit_id += total_sucessful_extention_for_subject;

        }
    }

    gpu_host_mem_free((void**)&local_h_GPUBlastInitHitList);
    local_h_GPUBlastInitHitList = GPUBlastInitHitList_tmp;
    GPUBlastInitHitList_tmp = NULL;

}

void autoset_successful_extension_in_volume(const Int4 num_sequences)
{
    int i;
    long int result = 0;
    for(i=0; i<num_sequences; i++)
    {
        result += local_h_Hits[i*4+1];
    }
    successful_extension_in_volume = result;

}

void H_BLAST_get_data( Int4* h_Hits,
                          const Int4 h_Hits_bytes,
                          GPUBlastInitHitList* h_GPUBlastInitHitList,
                          const Int4 h_GPUBlastInitHitList_bytes,
                          const Int4 num_sequences)

{
    int gpu_device_count=num_of_gpu_cards_used;

    for(int i=0; i<gpu_device_count; i++)
    {
        gpu_stream_synchronize(i);
        ( gpu_memcpy( &(local_h_Hits[stream_blocks[i].h_Hits_offset]), stream_blocks[i].d_Hits,
                      stream_blocks[i].h_Hits_bytes, cudaMemcpyDeviceToHost ) );
        ( gpu_memcpy( &(local_h_GPUBlastInitHitList[stream_blocks[i].h_GPUBlastInitHitList_offset]), stream_blocks[i].d_GPUBlastInitHitList,
                      stream_blocks[i].h_GPUBlastInitHitList_bytes, cudaMemcpyDeviceToHost));
    }

    set_max_gpu_walltime();

    memory_defrag(h_Hits_bytes, h_GPUBlastInitHitList_bytes, num_sequences);

    autoset_successful_extension_in_volume(num_sequences);

    memcpy(  h_Hits, local_h_Hits, h_Hits_bytes ) ;
    memcpy( h_GPUBlastInitHitList, local_h_GPUBlastInitHitList, h_GPUBlastInitHitList_bytes ) ;

    if(local_h_Hits!=NULL) gpu_host_mem_free((void**)&local_h_Hits);
    free(local_h_GPUBlastInitHitList);
    local_h_GPUBlastInitHitList=NULL;
}

void H_BLAST_free_memory()
{

    int gpu_device_count=num_of_gpu_cards_used;


    for(int cuda_device_id=0; cuda_device_id<gpu_device_count; cuda_device_id++)
    {
        gpu_set_device(stream_blocks[cuda_device_id].GPU_device_id);

        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_Hits ) ;
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_GPUBlastInitHitList ) ;
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_Database2Dpadded );
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_RepeatedSubstitutionMatrix );
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_RepeatedPV );
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_ThickBackbone );
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_overflow );
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_RepeatedDiag_array );
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_hash_table );
        gpu_device_mem_free( (void**)&stream_blocks[cuda_device_id].d_one_hit_buffer );

    }

}


Boolean H_BLAST_check_memory(const LookupTableWrap* lookup_wrap,
                                const Blast_ExtendWord* ewp,
                                const BlastGPUOptions* gpu_options,
                                const int h_Database2Dpadded_bytes,
                                const Int4 h_Hits_bytes, const Int4 h_GPUBlastInitHitList_bytes,
                                const Int4 Group_number, const Int4 num_queries,
                                const Int4 query_length,
                                const char mem_ecc_opton
                               )
{

    if( num_queries > NUM_QUERIES_MAX )
        return FALSE;
    if( query_length > QUERY_LENGTH_MAX )
        return FALSE;

    int h_RepeatedSubstitutionMatrix_bytes = SUBSTITUTION_MATRIX_LENGTH * (gpu_options->num_blocksx) * sizeof(char);

    BlastAaLookupTable *lookup = (BlastAaLookupTable*) lookup_wrap->lut;

    int pv_length = (lookup->backbone_size >> PV_ARRAY_BTS) + 1;
    Round2Multiple(&pv_length);
    int h_RepeatedPV_bytes = (gpu_options->num_blocksx) * pv_length * sizeof(PV_ARRAY_TYPE);

    int h_ThickBackbone_bytes = (lookup->backbone_size +1) * sizeof(AaLookupSmallboneCell);

    int h_overflow_bytes = (lookup->overflow_size+1) * sizeof(Uint2);

    int diag_array_length = ewp->diag_table->diag_array_length;
    unsigned long long h_RepeatedDiag_array_bytes =
        (unsigned long long)diag_array_length * (gpu_options->num_blocksx)*(gpu_options->num_threadsx)*sizeof(short);

    if(0!=get_gpu_kernel_option())
    {
        h_RepeatedDiag_array_bytes /=4;
        h_RepeatedDiag_array_bytes += ONE_HIT_BUFFER_SIZE*(gpu_options->num_blocksx)*(gpu_options->num_threadsx) * sizeof(int)/4;
    }

    if(0==num_of_gpu_cards_used)
    {
        return FALSE;
    }

    int acutal_h_Database2Dpadded_bytes = h_Database2Dpadded_bytes / num_of_gpu_cards_used;
    if(num_of_gpu_cards_used>1)
        acutal_h_Database2Dpadded_bytes *=IMB_factor_of_gpudb_group;

    unsigned long long Total_global_memory_bytes  =
        acutal_h_Database2Dpadded_bytes + gpu_options->num_blocksx * gpu_options->num_threadsx * sizeof(int) +
        h_RepeatedSubstitutionMatrix_bytes +
        h_Hits_bytes +
        h_RepeatedPV_bytes +
        h_ThickBackbone_bytes +
        h_overflow_bytes +
        h_RepeatedDiag_array_bytes +
        h_GPUBlastInitHitList_bytes;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    unsigned long long GPU_global_memory = (1==mem_ecc_opton)?(unsigned long long) (deviceProp.totalGlobalMem*5.0/6.0):(unsigned long long) (deviceProp.totalGlobalMem);

    if( Total_global_memory_bytes > GPU_global_memory)
    {
        printf("Not enough mem on the GPU: %ld bytes ,just has %ld \n", Total_global_memory_bytes, GPU_global_memory);
        return FALSE;
    }

    return TRUE;
}

Boolean GPU_BLAST_check_availability()
{

    int deviceCount = 0;
    if (cudaGetDeviceCount(&deviceCount) != 0)
        return FALSE;

    if ( 0 == deviceCount )
        return FALSE;

    return TRUE;
}
