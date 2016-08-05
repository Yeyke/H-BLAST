#ifndef H_BLAST_KERNEL_H
#define H_BLAST_KERNEL_H

#define ONE_HIT_BUFFER_SIZE 512

typedef struct ExtendLeftReturn {

    int maxscore;
    int length;

} ExtendLeftReturn;


typedef struct ExtendRightReturn {

    int s_last_off;
    int maxscore;
    int length;

} ExtendRightReturn;


typedef struct ExtendTwoHitReturn {

    int s_last_off;
    int hsp_q;
    int hsp_s;
    int hsp_len;
    int MAX_score;
    int right_extend;

    int left_score;
    int left_disp;
    int right_disp;

} ExtendTwoHitReturn;


#define AA_HITS_PER_CELL 3 /**< maximum number of hits in one lookup table cell */

__global__ void
H_BLAST_kernel_TwoHit1(
                          const int *d_RepeatedSubstitutionMatrix,
                          const int SubstitutionMatrix_length,
                          const PV_ARRAY_TYPE* d_RepeatedPV,
                          const AaLookupSmallboneCell* d_ThickBackbone,
                          const unsigned short* d_overflow,
                          const int* d_RepeatedSequenceMaxlength_vector,
                          const int SequenceMaxlength_vector_stride,
                          const int NumSequences,
                          const int pv_length,
                          const int Query_length,
                          const int diag_mask,
                          const int diag_array_length,
                          const int num_queries,
                          const int word_length,
                          const int index_mask,
                          int diag_offset,
                          int* d_Hits,
                          unsigned int total_hit_list_length);


__global__ void
H_BLAST_kernel_TwoHit2(
                          const int *d_RepeatedSubstitutionMatrix,
                          const int SubstitutionMatrix_length,
                          const PV_ARRAY_TYPE* d_RepeatedPV,
                          const AaLookupSmallboneCell* d_ThickBackbone,
                          const unsigned short* d_overflow,
                          const int* d_RepeatedSequenceMaxlength_vector,
                          const int SequenceMaxlength_vector_stride,
                          const int NumSequences,
                          const int pv_length,
                          const int Query_length,
                          const int diag_mask,
                          const int diag_array_length,
                          const int num_queries,
                          const int word_length,
                          const int index_mask,
                          int diag_offset,
                          int* d_Hits,
                          unsigned int total_hit_list_length,
                          int* hash_table,
                          int* d_one_hit_buffer);


static void
__device__ s_BlastAaExtendTwoHit(GPUBlastInitHitList* data,
                                 const char* matrix,
                                 const int* subject,
                                 const int subject_length,
                                 const char* query,
                                 const int query_length,
                                 const int s_left_off,
                                 int s_right_off,
                                 int q_right_off,
                                 const int dropoff,
                                 const int word_size,
                                 int* s_last_off,
                                 int* right_extend,
                                 char* AA_shared
                                 );


__device__ ExtendLeftReturn s_BlastAaExtendLeft(const char* matrix,
                                                const int* subject,
                                                const char * query,
                                                const int s_off,
                                                const int q_off,
                                                const int dropoff,
                                                int maxscore);

__device__ ExtendRightReturn s_BlastAaExtendRight(const char* matrix,
                                                  const int* subject,
                                                  const char* query,
                                                  const int s_off,
                                                  const int q_off,
                                                  const int dropoff,
                                                  int maxscore,
                                                  const int Sequence_actual_length,
                                                  const int Query_length);

#endif /* H_BLAST_KERNEL_H */
