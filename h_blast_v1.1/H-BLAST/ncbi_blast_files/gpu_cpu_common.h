#ifndef GPU_CPU_COMMON_H
#define GPU_CPU_COMMON_H

#include <algo/blast/core/blast_aalookup.h>
#include <algo/blast/core/lookup_wrap.h>
#include <algo/blast/core/blast_extend.h>
#include <algo/blast/core/blast_nalookup.h>


#define NUM_SEQUENCES                   99
#define PERCENTAGE                      99
#define NUM_THREADSX                   192
#define NUM_THREADSX_MAX               65536
#define NUM_BLOCKSX                    512
#define NUM_BLOCKSX_MAX              65535
#define GPU                              1
#define GPU_VERBOSE                     0
#define METHOD                           1
#define PITCH                           16
#define CPU_THREADS                      1

#define QUERY_LENGTH_MAX             32768
#define NUM_QUERIES_MAX                500

#define SUBSTITUTION_MATRIX_LENGTH     896
#define H_HITS_SIZE_BYTE 4
#define BRIEF_INFO 1

const static float IMB_factor_of_gpudb_group = 1.1;

typedef struct GPUBlastInitHitList
{

    unsigned short hsp_q;
    unsigned short hsp_s;
    unsigned short offset_from_first_position;
    unsigned short hsp_len;
    unsigned short score;

} GPUBlastInitHitList;


#ifdef __cplusplus
extern "C" {
#endif


    void  H_BLAST_pre(
        const BLAST_SequenceBlk* query, const BlastQueryInfo* query_info,
        const LookupTableWrap* lookup_wrap,
        const Blast_ExtendWord* ewp, const BlastInitialWordParameters* word_params,
        const BlastGPUOptions* gpu_options,  Int4* Sequence_length_vector,
        const int h_Database2Dpadded_bytes,
        const Int4 Group_number,
        const Int4 h_Hits_bytes, const Int4 h_GPUBlastInitHitList_bytes);

    void H_BLAST_get_data( Int4* h_Hits,
                           const Int4 h_Hits_bytes,
                           GPUBlastInitHitList* h_GPUBlastInitHitList,
                           const Int4 h_GPUBlastInitHitList_bytes,
                           const Int4 num_sequences);


    void Read_H_BLAST_GPU_DB(const char* GPU_volume_name,
                             const unsigned long h_H_BLAST_database_bytes,
                             const int volume_id);

    void all_gpu_abnormally_exit(int exit_code);
    void all_gpu_quit();
    void gpu_quit();
    void gpu_abnormally_exit(int exit_code);

    double get_gpu_walltime();
    void set_max_gpu_walltime();
    void init_max_and_current_gpu_walltime();

    void set_gpu_kernel_option(int flag);

    void set_first_gpu_id(unsigned int id);


    void  H_BLAST(
        const int cuda_device_id,
        const BLAST_SequenceBlk* query, const BlastQueryInfo* query_info,
        const LookupTableWrap* lookup_wrap,
        const Blast_ExtendWord* ewp, const BlastInitialWordParameters* word_params,
        const BlastGPUOptions* gpu_options, const int* Sequence_length_vector
    );


    void H_BLAST_free_memory();

    Boolean H_BLAST_check_memory(const LookupTableWrap* lookup_wrap,
                                 const Blast_ExtendWord* ewp,
                                 const BlastGPUOptions* gpu_options,
                                 const int h_Database2Dpadded_bytes,
                                 const Int4 h_Hits_bytes, const Int4 h_GPUBlastInitHitList_bytes,
                                 const Int4 Group_number, const Int4 num_queries,
                                 const Int4 query_length,
                                 const char mem_ecc_opton
                                );

    Boolean GPU_BLAST_check_availability();
    void Round2Multiple(Int4* MaxLength);

    void SetNumOfGpuCards(const int expected_num_of_cards);

    int get_num_of_gpu_cards_used();

    void set_MAX_SUCCESSFUL_HITS_PER_SEQUENCE(const unsigned int val);
    unsigned int get_MAX_SUCCESSFUL_HITS_PER_SEQUENCE();

    void lut_sorting(BlastAaLookupTable *lookup);

    long int get_successful_extension_in_volume () ;

    void Clean_H_BLAST_database();

    int get_processor_no();
    int get_core_count_per_processor();



#ifdef __cplusplus
}
#endif

#endif /* GPU_CPU_COMMON_H */
