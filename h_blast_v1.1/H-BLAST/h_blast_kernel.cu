#include <algo/blast/core/gpu_cpu_common.h>
#include "h_blast_kernel.h"
#define CHARSIZE 5

extern __shared__ int array[];
__constant__ char d_Query_head_const[4];
__constant__ char d_Query_const[40960];
__constant__ short d_x_dropoff_const[900];
__constant__ short d_cutoff_score_const[900];
__constant__ short d_query_offset_const[900];
__constant__ int d_length_vector_const[256];

__constant__ int* d_Database2Dpadded_const;
__constant__ short* d_diag_array_const;
__constant__ GPUBlastInitHitList* d_GPUBlastInitHitList_const;
__constant__ int hit_list_base_const;
__constant__ int diag_offset_const;



__global__ void
H_BLAST_kernelTwoHit2(
    const int* d_RepeatedSubstitutionMatrix,
    const int SubstitutionMatrix_length,
    const PV_ARRAY_TYPE* d_RepeatedPV,
    const AaLookupSmallboneCell* d_ThickBackbone,
    const unsigned short* d_overflow,
    const int Sequence_length_vector_stride,
    const int NumSequences,
    const int pv_length,
    const int Query_length,
    const int diag_mask,
    const int diag_array_length,
    const int num_queries,
    const int word_length,
    const int index_mask,
    const int window,
    int* d_Hits,
    unsigned int total_hit_list_length,
    int* hash_table,
    int* d_one_hit_buffer)
{
    int* d_Database2Dpadded = d_Database2Dpadded_const;
    short* d_diag_array = d_diag_array_const;
    GPUBlastInitHitList* d_GPUBlastInitHitList = d_GPUBlastInitHitList_const;
    int diag_offset;


    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int blockDimx = blockDim.x;

    int* SubstitutionMatrix_shared_int = (int*)array;
    int* result_space_used = (int*)  &SubstitutionMatrix_shared_int[SubstitutionMatrix_length>>2];
    char* Sequence_AA = (char*) &result_space_used[blockDimx];
    int* hash_table_shared = (int*)( &Sequence_AA[5* blockDimx * 4] );
    PV_ARRAY_TYPE* pv_shared = (PV_ARRAY_TYPE*) (&hash_table_shared[16]);

    int pv_Index_global = bx*pv_length + tx;
    int pv_Index_shared = tx;
    while( pv_Index_shared < pv_length )
    {
        pv_shared[ pv_Index_shared ] = d_RepeatedPV[ pv_Index_global ];

        pv_Index_global += blockDimx;
        pv_Index_shared += blockDimx;
    }

    __syncthreads();

    int SubstitutionMatrix_Index_global = bx*(SubstitutionMatrix_length>>2) + tx;
    int SubstitutionMatrix_Index_shared = tx;
    while( SubstitutionMatrix_Index_shared < (SubstitutionMatrix_length>>2) )
    {
        SubstitutionMatrix_shared_int[ SubstitutionMatrix_Index_shared ] = d_RepeatedSubstitutionMatrix[ SubstitutionMatrix_Index_global ];

        SubstitutionMatrix_Index_global += blockDimx;
        SubstitutionMatrix_Index_shared += blockDimx;
    }

    __syncthreads();

    int hash_id = tx;

    while(hash_id<16)
    {
        hash_table_shared[hash_id] = hash_table[hash_id];
        hash_id += blockDimx;
    }


    __syncthreads();

    char* SubstitutionMatrix_shared_char = (char*)SubstitutionMatrix_shared_int;

    int loop_counter = 0;

    d_one_hit_buffer = & (d_one_hit_buffer[ (bx*blockDimx*ONE_HIT_BUFFER_SIZE)/4 + (threadIdx.x/4)*ONE_HIT_BUFFER_SIZE]);

    Sequence_AA = &Sequence_AA[(threadIdx.x/4)*80];

    int total_threads = blockDim.x*gridDim.x;
    int left_room = total_hit_list_length/total_threads;

    int current_subject_offset = left_room*( bx*blockDimx + (tx &0xfe0))+ left_room * 8* ((tx/4)&0x3);

    left_room *= 8;

    int hit_list_base = hit_list_base_const;

    if( 0==(threadIdx.x&0x10) )
    {
        d_GPUBlastInitHitList = &(d_GPUBlastInitHitList_const[current_subject_offset]);
        hit_list_base += current_subject_offset;
    }
    else
    {
        d_GPUBlastInitHitList = &(d_GPUBlastInitHitList_const[current_subject_offset+left_room-1]);
        hit_list_base += current_subject_offset+left_room-1;
    }


    for(loop_counter=0; loop_counter<4; loop_counter++)
    {
        int init_subject_id_in_block = (threadIdx.x /16) *16 + (threadIdx.x&0xf)/4 + 4*loop_counter;

        int MBound = (total_hit_list_length/(blockDim.x*gridDim.x))*2*(3-loop_counter);

        int Sequence = (bx * blockDimx + init_subject_id_in_block);

        GPUBlastInitHitList data;

        int Group_size = gridDim.x*blockDimx;

        if(0==loop_counter)
        {
            int diag_array_base = ((bx*blockDimx/4+threadIdx.x/4)*diag_array_length);
            d_diag_array = &d_diag_array_const[diag_array_base];
        }
        else
        {
            int* d_diag_array_int = (int*) d_diag_array;

            for(int i= diag_array_length/2 - 4 + (threadIdx.x&0x3); i>=0; i-=4 )
                d_diag_array_int[i] = 0;
        }
        diag_offset = diag_offset_const;

        d_Database2Dpadded = &(d_Database2Dpadded_const[init_subject_id_in_block+ blockDim.x]);

        int S_AA_buffer;

        while( Sequence < NumSequences )
        {
            int Hits = 0;
            int s_last_off = 0;
            int hits_extended = 0;
            int successful_extensions = 0;
            int Sequence_length = d_length_vector_const[ Sequence/Group_size ];

            int block_base = (bx*blockDimx*Sequence_length)>>2;
            d_Database2Dpadded = & (d_Database2Dpadded[ block_base ]);
            ((int*) Sequence_AA)[0] = d_Database2Dpadded[ 0 ];
            S_AA_buffer = d_Database2Dpadded[ blockDimx ];

            result_space_used[threadIdx.x] = 0;
            ((int*)Sequence_AA)[2] = 0;

            int index = (Sequence_AA[ 0 ] << CHARSIZE) |  Sequence_AA[ 1 ];
            int update_offset_index = 0;

            for(int s_off = 0; (s_off <= Sequence_length - word_length); ++s_off)
            {
                if(0==index)
                    break;

                if ( ((s_off&0x3)==0) && (s_off < Sequence_length-4) )
                {

                    if (0==(threadIdx.x & 0x3))
                    {
                        ((int*)Sequence_AA)[((s_off>>2)+1)&0x1] = S_AA_buffer;
                        S_AA_buffer = d_Database2Dpadded[ (((s_off>>2) + 2)*blockDimx) ];

                    }

                }

                index = ((index<<CHARSIZE) | Sequence_AA[ (s_off+2) &0x7 ]) & index_mask;

                int numhits = 0;


                unsigned short to_be_extension_count = 0;

                int query_offset, diag_coord, query_offset_buffer1, query_offset_buffer2;
                short diagonal, diagonal_buffer;
                short last_hit;
                short flag;

                if( PV_TEST(pv_shared, index, PV_ARRAY_BTS) )
                {
                    numhits = d_ThickBackbone[index].num_used;
                    Hits += numhits;

                    const unsigned short *src;
                    if( numhits <= AA_HITS_PER_CELL )
                        src = d_ThickBackbone[index].payload.entries;
                    else
                        src = &(d_overflow[d_ThickBackbone[index].payload.overflow_cursor]);

                    query_offset_buffer1 = src[(threadIdx.x&0x3)];
                    query_offset_buffer2 = src[4+(threadIdx.x&0x3)];
                    numhits -=4;
                    diag_coord = (query_offset_buffer1 - s_off) & diag_mask;
                    diagonal_buffer = d_diag_array[diag_coord];

                    int i;
                    for(i = (threadIdx.x&0x3); i < numhits; i+=4)
                    {

                        query_offset = query_offset_buffer1;
                        diagonal = diagonal_buffer;

                        diagonal_buffer = d_diag_array[(query_offset_buffer2 - s_off) & diag_mask];
                        query_offset_buffer1= query_offset_buffer2;
                        query_offset_buffer2 = src[i+8];

                        diag_coord = (query_offset - s_off) & diag_mask;
                        last_hit = diagonal & 0x7FFF;
                        flag = diagonal & 0x8000;


                        if (flag)
                        {
                            if ( (s_off + diag_offset) >= last_hit)
                                d_diag_array[diag_coord] = (s_off + diag_offset) & 0x7FFF;

                        }
                        else
                        {
                            last_hit = last_hit - diag_offset;
                            int diff = s_off - last_hit;

                            if (diff >= window)
                            {
                                d_diag_array[diag_coord] = (s_off + diag_offset) & 0x7FFF;
                            }
                            else if (diff >= word_length)
                            {
                                successful_extensions = ((int*)Sequence_AA)[2];
                                to_be_extension_count = (unsigned short) (successful_extensions & 0xffff);
                                successful_extensions >>= 16;

                                d_diag_array[diag_coord] = (s_off + diag_offset) & 0x7FFF;
                                int thread_interaction_signal = (__ballot(1)>>(threadIdx.x&0x1c))&0xf;
                                int tmpl = hash_table_shared[thread_interaction_signal];
                                int extension_invoked = ( tmpl >> 9)&0x7 ;
                                tmpl = (tmpl>>((threadIdx.x&0x3)*3 ))&0x7;

                                if(to_be_extension_count + extension_invoked <= ONE_HIT_BUFFER_SIZE)
                                {
                                    if (to_be_extension_count + tmpl <= 17 )
                                    {
                                        ((int*)Sequence_AA)[2+to_be_extension_count + tmpl] = query_offset + last_hit*65536;
                                    }
                                    else
                                    {
                                        d_one_hit_buffer[to_be_extension_count + tmpl-18] = query_offset + last_hit*65536;
                                    }
                                    to_be_extension_count+= extension_invoked;
                                    ((int*)Sequence_AA)[2] = successful_extensions*65536 + to_be_extension_count;
                                }

                            }

                        }
                    }

                    if(i<(numhits+4))
                    {
                        diagonal = diagonal_buffer;
                        query_offset = query_offset_buffer1;
                        diag_coord = (query_offset - s_off) & diag_mask;
                        last_hit = diagonal & 0x7FFF;
                        flag = diagonal & 0x8000;

                        if (flag)
                        {
                            if ( (s_off + diag_offset) >= last_hit)
                                d_diag_array[diag_coord] = (s_off + diag_offset) & 0x7FFF;

                        }
                        else
                        {
                            last_hit = last_hit - diag_offset;
                            int diff = s_off - last_hit;

                            if (diff >= window)
                            {
                                d_diag_array[diag_coord] = (s_off + diag_offset) & 0x7FFF;
                            }

                            else if (diff >= word_length)
                            {
                                successful_extensions = ((int*)Sequence_AA)[2];
                                to_be_extension_count = (unsigned short) (successful_extensions & 0xffff);
                                successful_extensions >>= 16;

                                d_diag_array[diag_coord] = (s_off + diag_offset) & 0x7FFF;


                                int thread_interaction_signal = (__ballot(1)>>(threadIdx.x&0x1c))&0xf;
                                int tmpl = hash_table_shared[thread_interaction_signal];
                                int extension_invoked = ( tmpl >> 9)&0x7 ;
                                tmpl = (tmpl>>((threadIdx.x&0x3)*3 ))&0x7;

                                if(to_be_extension_count + extension_invoked <= ONE_HIT_BUFFER_SIZE)
                                {
                                    if (to_be_extension_count + tmpl <= 17 )
                                    {
                                        ((int*)Sequence_AA)[2+to_be_extension_count + tmpl] = query_offset + last_hit*65536;
                                    }
                                    else
                                    {
                                        d_one_hit_buffer[to_be_extension_count + tmpl-18] = query_offset + last_hit*65536;
                                    }
                                    to_be_extension_count+= extension_invoked;
                                    ((int*)Sequence_AA)[2] = successful_extensions*65536 + to_be_extension_count;

                                }

                            }
                        }
                    }
                }

                successful_extensions = ((int*)Sequence_AA)[2];
                to_be_extension_count = (unsigned short) (successful_extensions & 0xffff);
                successful_extensions >>= 16;

                if( (1==(s_off&0x1)) && (to_be_extension_count>0) )
                {
                    int s_off_previous = s_off-1;
                    for(int i=(threadIdx.x&0x3); i<to_be_extension_count; i+=4)
                    {

                        if(i>=update_offset_index)
                            s_off_previous=s_off;

                        if(i<17)
                        {
                            query_offset = ((unsigned short*)Sequence_AA)[6+2*i];
                            last_hit = ((unsigned short*)Sequence_AA)[7+2*i];
                        }
                        else
                        {
                            int tmp = d_one_hit_buffer[i-17];
                            query_offset = tmp&0xffff;
                            last_hit = tmp>>16;
                        }

                        diag_coord = (query_offset - s_off_previous) & diag_mask;

                        successful_extensions = ((int*)Sequence_AA)[2];
                        successful_extensions >>= 16;

                        int left = 0;
                        int right = num_queries;
                        int mid=0;

                        while((right -left) > 1)
                        {
                            mid = (left+right)>>1;
                            if ( ((int)d_query_offset_const[mid]) <= query_offset )
                                left = mid;
                            else
                                right = mid;
                        }
                        int k = left;

                        int right_extend;
                        s_BlastAaExtendTwoHit(&data,
                                              SubstitutionMatrix_shared_char,
                                              (d_Database2Dpadded),
                                              Sequence_length,
                                              d_Query_const,
                                              Query_length,
                                              last_hit + word_length,
                                              s_off_previous,
                                              query_offset,
                                              (int)d_x_dropoff_const[k],
                                              word_length,
                                              &s_last_off,&right_extend,
                                              Sequence_AA
                                             );

                        ++hits_extended;

                        if ( data.score >= (int)d_cutoff_score_const[k] )
                        {
                            data.offset_from_first_position = s_off_previous - data.hsp_s;

                            int thread_interaction_signal = (__ballot(1)>>(threadIdx.x&0x1c))&0xf;
                            int tmpl = hash_table_shared[thread_interaction_signal];
                            int storage_invoked_count = ( tmpl >> 9)&0x7 ;
                            tmpl = ((tmpl>>((threadIdx.x&0x3)*3 ))&0x7) -1;

                            if(left_room >= (successful_extensions + storage_invoked_count ) )
                            {
                                if( 0==(threadIdx.x&0x10) )
                                {
                                    d_GPUBlastInitHitList[successful_extensions+tmpl] = data;
                                }
                                else
                                {
                                    d_GPUBlastInitHitList[0-successful_extensions-tmpl] = data;
                                }

                            }
                            successful_extensions+= storage_invoked_count ;
                            ((int*)Sequence_AA)[2] = successful_extensions*65536 + to_be_extension_count;

                        }

                        if (right_extend)
                        {

                            d_diag_array[diag_coord] = ( (s_last_off - (word_length - 1) + diag_offset) & 0x7FFF ) | 0x8000;

                        }

                        result_space_used[threadIdx.x] = 0;
                    }

                    if( (0==(threadIdx.x&0x3))  )
                    {
                        ((int*)Sequence_AA)[2] = ((int*)Sequence_AA)[2] & 0xffff0000;
                        update_offset_index = 0;
                    }
                }
                else
                {
                    update_offset_index = to_be_extension_count;
                }
            }


            result_space_used[threadIdx.x] = hits_extended;
            hits_extended += result_space_used[threadIdx.x^0x1];

            result_space_used[threadIdx.x] = hits_extended;
            hits_extended += result_space_used[threadIdx.x^0x2];

            successful_extensions = ((int*)Sequence_AA)[2];
            successful_extensions >>= 16;


            diag_offset += Sequence_length + window;


            int iTmp = Sequence*4;

            result_space_used[threadIdx.x] = successful_extensions;

            int iTmp2 = result_space_used[threadIdx.x ^ 0x10];

            if ( ( successful_extensions > left_room ) || ((left_room - MBound ) < (successful_extensions + iTmp2)) )
            {
                d_Hits[ iTmp+3 ] = -1;
            }
            else
            {
                d_Hits[ iTmp+3 ] = hit_list_base;

                if(threadIdx.x & 0x10)
                {
                    d_GPUBlastInitHitList = & (d_GPUBlastInitHitList[0-successful_extensions]);
                    hit_list_base -= successful_extensions;
                }
                else
                {
                    d_GPUBlastInitHitList = & (d_GPUBlastInitHitList[successful_extensions]);
                    hit_list_base += successful_extensions;
                }

                left_room -= iTmp2 + successful_extensions;
            }

            d_Hits[ iTmp ] = Hits;
            d_Hits[ iTmp + 1 ] = successful_extensions;
            d_Hits[ iTmp + 2 ] = hits_extended;

            iTmp = (Group_size - bx*blockDimx) * (Sequence_length>>2);

            Sequence += Group_size;

            d_Database2Dpadded = &(d_Database2Dpadded[iTmp]);

        }
    }
}

__global__ void
H_BLAST_kernelTwoHit1(
    const int* d_RepeatedSubstitutionMatrix,
    const int SubstitutionMatrix_length,
    const PV_ARRAY_TYPE* d_RepeatedPV,
    const AaLookupSmallboneCell* d_ThickBackbone,
    const unsigned short* d_overflow,
    const int Sequence_length_vector_stride,
    const int NumSequences,
    const int pv_length,
    const int Query_length,
    const int diag_mask,
    const int diag_array_length,
    const int num_queries,
    const int word_length,
    const int index_mask,
    const int window,
    int* d_Hits,
    unsigned int total_hit_list_length)
{
    int* d_Database2Dpadded = d_Database2Dpadded_const;
    short* d_diag_array = d_diag_array_const;
    GPUBlastInitHitList* d_GPUBlastInitHitList = d_GPUBlastInitHitList_const;
    int hit_list_base = hit_list_base_const;
    int diag_offset = diag_offset_const;

    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int blockDimx = blockDim.x;

    int* SubstitutionMatrix_shared_int = (int*)array;
    int* result_space_used = (int*)  &SubstitutionMatrix_shared_int[SubstitutionMatrix_length>>2];
    char* Sequence_AA = (char*) &result_space_used[blockDimx];

    PV_ARRAY_TYPE* pv_shared = (PV_ARRAY_TYPE*) &Sequence_AA[5* blockDimx * 4];

    int pv_Index_global = bx*pv_length + tx;
    int pv_Index_shared = tx;
    while( pv_Index_shared < pv_length )
    {
        pv_shared[ pv_Index_shared ] = d_RepeatedPV[ pv_Index_global ];

        pv_Index_global += blockDimx;
        pv_Index_shared += blockDimx;
    }

    __syncthreads();

    int SubstitutionMatrix_Index_global = bx*(SubstitutionMatrix_length>>2) + tx;
    int SubstitutionMatrix_Index_shared = tx;
    while( SubstitutionMatrix_Index_shared < (SubstitutionMatrix_length>>2) )
    {
        SubstitutionMatrix_shared_int[ SubstitutionMatrix_Index_shared ] = d_RepeatedSubstitutionMatrix[ SubstitutionMatrix_Index_global ];

        SubstitutionMatrix_Index_global += blockDimx;
        SubstitutionMatrix_Index_shared += blockDimx;
    }

    __syncthreads();

    char* SubstitutionMatrix_shared_char = (char*)SubstitutionMatrix_shared_int;

    {
        int total_threads = blockDim.x*gridDim.x;
        int left_room = total_hit_list_length/total_threads;
        int current_subject_offset = left_room*( bx*blockDimx + (tx&0xfe0))+ left_room * 2* (tx&0xf);
        left_room *= 2;

        if( 0==(tx&0x10) )
        {
            d_GPUBlastInitHitList = &(d_GPUBlastInitHitList[current_subject_offset]);
            hit_list_base += current_subject_offset;
        }
        else
        {
            d_GPUBlastInitHitList = &(d_GPUBlastInitHitList[current_subject_offset+left_room-1]);
            hit_list_base += current_subject_offset+left_room-1;
        }

        int Sequence = (bx * blockDimx + tx);

        GPUBlastInitHitList data;

        int Group_size = gridDim.x*blockDimx;
        int diag_array_base = ((bx*blockDimx+tx)*diag_array_length);

        Sequence_AA = &Sequence_AA[tx*20];
        d_diag_array = &d_diag_array[diag_array_base];

        d_Database2Dpadded = &(d_Database2Dpadded[tx+ blockDim.x]);

        int S_AA_buffer;

        while( Sequence < NumSequences )
        {
            int Hits = 0;
            int s_last_off = 0;
            int hits_extended = 0;
            int successful_extensions = 0;
            int Sequence_length = d_length_vector_const[ Sequence/Group_size ];

            int block_base = (bx*blockDimx*Sequence_length)>>2;
            d_Database2Dpadded = & (d_Database2Dpadded[ block_base ]);
            ((int*) Sequence_AA)[0] = d_Database2Dpadded[ 0 ];
            S_AA_buffer = d_Database2Dpadded[ blockDimx ];

            int index = (Sequence_AA[ 0 ] << CHARSIZE) |  Sequence_AA[ 1 ];

            int numhits_buffer = 0;
            const unsigned short *src_buffer=NULL;
            index = ((index<<CHARSIZE) | Sequence_AA[ 2 ]) & index_mask;
            short s_off_mask = diag_offset;
            if( PV_TEST(pv_shared, index, PV_ARRAY_BTS) )
            {
                numhits_buffer = d_ThickBackbone[index].num_used;

                if( numhits_buffer <= AA_HITS_PER_CELL )
                    src_buffer = d_ThickBackbone[index].payload.entries;
                else
                    src_buffer = &(d_overflow[d_ThickBackbone[index].payload.overflow_cursor]);
            }

            int to_be_extension_count = 0;
            int update_offset_index = 0;
            int s_off, query_offset,  diag_coord[2];
            short last_hit;

            for(s_off = 0; (s_off <= Sequence_length - word_length); ++s_off, s_off_mask++)
            {

                if ( ((s_off&0x3)==0) && (s_off < Sequence_length-4) )
                {
                    ((int*)Sequence_AA)[((s_off>>2)+1)&0x1] = S_AA_buffer;
                    S_AA_buffer = d_Database2Dpadded[ (((s_off>>2) + 2)*blockDimx) ];

                    if(0==index)
                        break;
                }

                int numhits = numhits_buffer;
                const unsigned short *src = src_buffer;

                if(s_off < (Sequence_length-word_length))
                {

                    index = ((index<<CHARSIZE) | Sequence_AA[ (s_off+3) &0x7 ]) & index_mask;

                    if( PV_TEST(pv_shared, index, PV_ARRAY_BTS) )
                    {
                        numhits_buffer = d_ThickBackbone[index].num_used;

                        if( numhits_buffer <= AA_HITS_PER_CELL )
                            src_buffer = d_ThickBackbone[index].payload.entries;
                        else
                            src_buffer = &(d_overflow[d_ThickBackbone[index].payload.overflow_cursor]);
                    }
                    else
                    {
                        numhits_buffer = 0;
                        src_buffer = NULL;
                    }


                }

                int query_offset_buffer1, query_offset_buffer2;
                short diagonal, diagonal_buffer;

                if(numhits>0)
                {

                    int diff;
                    query_offset_buffer1 = src[0];
                    query_offset_buffer2 = src[1];
                    diag_coord[0] = (query_offset_buffer1 - s_off) & diag_mask;
                    diag_coord[1] = (query_offset_buffer2 - s_off) & diag_mask;
                    Hits += numhits;
                    numhits --;

                    diagonal_buffer = d_diag_array[diag_coord[0]];


                    for(int i = 0; i < numhits; ++i)
                    {
                        query_offset = query_offset_buffer1;
                        diagonal = diagonal_buffer;
                        diagonal_buffer = d_diag_array[diag_coord[1]];
                        query_offset_buffer1= query_offset_buffer2;
                        last_hit = diagonal & 0x7FFF;
                        query_offset_buffer2 = src[i+2];

                        diff = (s_off_mask  - last_hit) | (diagonal & 0x8000);
                        last_hit = last_hit - diag_offset;

                        if (diff>=word_length)
                        {
                            d_diag_array[diag_coord[0]] = s_off_mask;

                            if (diff < window)
                            {

                                if(to_be_extension_count==3)
                                {
                                    int left = 0;
                                    int right = num_queries;
                                    int mid=0;

                                    while((right -left) > 1)
                                    {
                                        mid = (left+right)>>1;
                                        if (((int)d_query_offset_const[mid]) <= query_offset )
                                            left = mid;
                                        else
                                            right = mid;
                                    }
                                    int k = left;

                                    int right_extend;
                                    s_BlastAaExtendTwoHit(&data,
                                                          SubstitutionMatrix_shared_char,
                                                          (d_Database2Dpadded),
                                                          Sequence_length,
                                                          d_Query_const,
                                                          Query_length,
                                                          last_hit + word_length,
                                                          s_off,
                                                          query_offset,
                                                          (int)d_x_dropoff_const[k],
                                                          word_length,
                                                          &s_last_off,
                                                          &right_extend,
                                                          Sequence_AA
                                                         );

                                    ++hits_extended;

                                    if ( data.score >= (int)d_cutoff_score_const[k] )
                                    {
                                        data.offset_from_first_position = s_off - data.hsp_s;

                                        if(left_room > successful_extensions)
                                        {
                                            if( 0==(threadIdx.x&0x10) )
                                                d_GPUBlastInitHitList[successful_extensions] = data;
                                            else
                                                d_GPUBlastInitHitList[0-successful_extensions] = data;
                                        }

                                        successful_extensions++;


                                    }


                                    if (right_extend)
                                    {

                                        d_diag_array[diag_coord[0]] = ( (s_last_off - (word_length - 1) + diag_offset) & 0x7FFF ) | 0x8000;

                                    }

                                }
                                else
                                {
                                    ((int*)Sequence_AA)[2+to_be_extension_count] = query_offset + last_hit*65536;
                                    to_be_extension_count++;
                                    if(0==(s_off&0x1))
                                    {
                                        update_offset_index = to_be_extension_count;
                                    }


                                }

                            }

                        }

                        diag_coord[0] = diag_coord[1];
                        diag_coord[1] = (query_offset_buffer2 - s_off) & diag_mask;
                    }

                    query_offset = query_offset_buffer1;
                    diagonal = diagonal_buffer;


                    last_hit = diagonal & 0x7FFF;
                    diff = (s_off_mask  - last_hit) | (diagonal & 0x8000);
                    last_hit = last_hit - diag_offset;


                    if (diff>=word_length)
                    {
                        d_diag_array[diag_coord[0]] = s_off_mask;

                        if (diff <window)
                        {
                            if(to_be_extension_count==3)
                            {
                                int left = 0;
                                int right = num_queries;
                                int mid=0;

                                while((right -left) > 1)
                                {
                                    mid = (left+right)>>1;
                                    if (((int)d_query_offset_const[mid]) <= query_offset )
                                        left = mid;
                                    else
                                        right = mid;
                                }
                                int k = left;


                                int right_extend;
                                s_BlastAaExtendTwoHit(&data,
                                                      SubstitutionMatrix_shared_char,
                                                      (d_Database2Dpadded),
                                                      Sequence_length,
                                                      d_Query_const,
                                                      Query_length,
                                                      last_hit + word_length,
                                                      s_off,
                                                      query_offset,
                                                      (int)d_x_dropoff_const[k],
                                                      word_length,
                                                      &s_last_off,
                                                      &right_extend,
                                                      Sequence_AA
                                                     );

                                ++hits_extended;

                                if ( data.score >= (int)d_cutoff_score_const[k] )
                                {
                                    data.offset_from_first_position = s_off - data.hsp_s;

                                    if(left_room > successful_extensions)
                                    {
                                        if( 0==(threadIdx.x&0x10) )
                                            d_GPUBlastInitHitList[successful_extensions] = data;
                                        else
                                            d_GPUBlastInitHitList[0-successful_extensions] = data;

                                    }
                                    successful_extensions++;
                                }

                                if (right_extend)
                                {
                                    d_diag_array[diag_coord[0]] = ( (s_last_off - (word_length - 1) + diag_offset) & 0x7FFF ) | 0x8000;
                                }

                            }
                            else
                            {
                                ((int*)Sequence_AA)[2+to_be_extension_count] = query_offset + last_hit*65536;
                                to_be_extension_count++;
                                if(0==(s_off&0x1))
                                {
                                    update_offset_index = to_be_extension_count;

                                }
                            }

                        }
                    }
                }

                if(1==(s_off&0x1))
                {
                    int s_off_previous = s_off-1;
                    for(int i=0; i<to_be_extension_count; i++)
                    {
                        if(i==update_offset_index)
                        {
                            s_off_previous ++;
                        }
                        query_offset = ((unsigned short*)Sequence_AA)[4+2*i];
                        last_hit = ((unsigned short*)Sequence_AA)[5+2*i];
                        diag_coord[0] = (query_offset - s_off_previous) & diag_mask;

                        int left = 0;
                        int right = num_queries;
                        int mid=0;

                        while((right -left) > 1)
                        {
                            mid = (left+right)>>1;
                            if (((int)d_query_offset_const[mid]) <= query_offset )
                                left = mid;
                            else
                                right = mid;
                        }
                        int k = left;


                        int right_extend;
                        s_BlastAaExtendTwoHit(&data,
                                              SubstitutionMatrix_shared_char,
                                              (d_Database2Dpadded),
                                              Sequence_length,
                                              d_Query_const,
                                              Query_length,
                                              last_hit + word_length,
                                              s_off_previous,
                                              query_offset,
                                              (int)d_x_dropoff_const[k],
                                              word_length,
                                              &s_last_off,&right_extend, Sequence_AA
                                             );

                        ++hits_extended;

                        if ( data.score >= (int)d_cutoff_score_const[k] )
                        {
                            data.offset_from_first_position = s_off_previous - data.hsp_s;

                            if(left_room > successful_extensions)
                            {
                                if( 0==(threadIdx.x&0x10) )
                                    d_GPUBlastInitHitList[successful_extensions] = data;
                                else
                                    d_GPUBlastInitHitList[0-successful_extensions] = data;
                            }

                            successful_extensions++;
                        }

                        if (right_extend)
                        {

                            d_diag_array[diag_coord[0]] = ( (s_last_off - (word_length - 1) + diag_offset) & 0x7FFF ) | 0x8000;
                        }
                    }

                    update_offset_index = to_be_extension_count = 0;

                }
            }

            if(to_be_extension_count>0)
            {
                int s_off_previous = s_off-1;
                for(int i=0; i<to_be_extension_count; i++)
                {
                    if(i==update_offset_index)
                    {
                        s_off_previous ++;
                    }
                    query_offset = ((unsigned short*)Sequence_AA)[4+2*i];
                    last_hit = ((unsigned short*)Sequence_AA)[5+2*i];
                    diag_coord[0] = (query_offset - s_off_previous) & diag_mask;

                    int left = 0;
                    int right = num_queries;
                    int mid=0;

                    while((right -left) > 1)
                    {
                        mid = (left+right)>>1;
                        if (((int)d_query_offset_const[mid]) <= query_offset )
                            left = mid;
                        else
                            right = mid;
                    }
                    int k = left;


                    int right_extend;
                    s_BlastAaExtendTwoHit(&data,
                                          SubstitutionMatrix_shared_char,
                                          (d_Database2Dpadded),
                                          Sequence_length,
                                          d_Query_const,
                                          Query_length,
                                          last_hit + word_length,
                                          s_off_previous,
                                          query_offset,
                                          (int)d_x_dropoff_const[k],
                                          word_length,
                                          &s_last_off,
                                          &right_extend,
                                          Sequence_AA
                                         );

                    ++hits_extended;

                    if ( data.score >= (int)d_cutoff_score_const[k] )
                    {
                        data.offset_from_first_position = s_off_previous - data.hsp_s;

                        if(left_room > successful_extensions)
                        {
                            if( 0==(threadIdx.x&0x10) )
                                d_GPUBlastInitHitList[successful_extensions] = data;
                            else
                                d_GPUBlastInitHitList[0-successful_extensions] = data;
                        }

                        successful_extensions++;


                    }

                    if (right_extend)
                    {

                        d_diag_array[diag_coord[0]] = ( (s_last_off - (word_length - 1) + diag_offset) & 0x7FFF ) | 0x8000;
                    }

                }

                update_offset_index = to_be_extension_count = 0;

            }



            diag_offset += Sequence_length + window;


            int iTmp = Sequence*4;

            result_space_used[threadIdx.x] = successful_extensions;

            int iTmp2 = result_space_used[threadIdx.x^0x10];

            if ( ( successful_extensions > left_room ) || ((left_room - successful_extensions) <  iTmp2) )
            {
                d_Hits[ iTmp+3 ] = -1;
            }
            else
            {
                d_Hits[ iTmp+3 ] = hit_list_base;

                if(threadIdx.x & 0x10)
                {
                    d_GPUBlastInitHitList = & (d_GPUBlastInitHitList[0-successful_extensions]);
                    hit_list_base -= successful_extensions;
                }
                else
                {
                    d_GPUBlastInitHitList = & (d_GPUBlastInitHitList[successful_extensions]);
                    hit_list_base += successful_extensions;
                }

                left_room -= iTmp2 + successful_extensions;
            }

            d_Hits[ iTmp ] = Hits;
            d_Hits[ iTmp + 1 ] = successful_extensions;
            d_Hits[ iTmp + 2 ] = hits_extended;

            iTmp = (Group_size - bx*blockDimx) * (Sequence_length>>2);

            Sequence += Group_size;

            d_Database2Dpadded = &(d_Database2Dpadded[iTmp]);

        }
    }
}


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
                                )
{

    int left_d = 0, right_d = 0;
    int left_score = 0, right_score = 0;
    int i, score = 0;


    for (i = 0; i < word_size; ++i)
    {
        score += matrix[ query[q_right_off+i] + (AA_shared[(s_right_off+i) &0x7]<<5) ];

        if (score > left_score)
        {
            left_score = score;
            right_d = i + 1;
        }
    }

    q_right_off += right_d;
    s_right_off += right_d;

    right_d = 0;
    (*right_extend) = FALSE;
    (*s_last_off) = s_right_off;

    ExtendLeftReturn ReturnLeftValue = s_BlastAaExtendLeft(matrix, subject, query, s_right_off-1, q_right_off-1, dropoff, 0);
    left_score = ReturnLeftValue.maxscore;
    left_d = ReturnLeftValue.length;

    ExtendRightReturn ReturnRightValue;
    if (left_d >= (s_right_off - s_left_off))
    {
        (*right_extend) = TRUE;

        ReturnRightValue = s_BlastAaExtendRight(matrix, subject, query, s_right_off, q_right_off, dropoff, left_score, subject_length, query_length);
        right_score = ReturnRightValue.maxscore;
        right_d = ReturnRightValue.length;
        (*s_last_off) = ReturnRightValue.s_last_off;
    }

    (*data).hsp_q = q_right_off - left_d;
    (*data).hsp_s = s_right_off - left_d;;
    (*data).hsp_len = left_d + right_d;
    (*data).score = MAX(left_score, right_score);

}


__device__ ExtendRightReturn s_BlastAaExtendRight(const char* matrix,
        const int* subject,
        const char* query,
        const int s_off,
        const int q_off,
        const int dropoff,
        int maxscore,
        const int Sequence_length,
        const int Query_length)
{
    int n, best_i = -1;
    int score = maxscore;

    int subject_buffer1, subject_buffer2;
    int temp = blockDim.x;
    int index = ((s_off)>>2)*temp;
    char query_temp_buffer = query[q_off];
    subject_buffer1 = subject[index];
    index += temp;
    subject_buffer2 = subject[index ];
    char* subject_buffer1_ptr = (char*) &subject_buffer1;


    n = MIN(Sequence_length - s_off, Query_length - q_off);

    char query_temp = 0, sequence_temp = 0,
                                         flag = 0;


    int temp2 = 4-(s_off&0x3);
    int n1=n-temp2;

    int n2 = n1>>2;
    int k=0;

    if(n>0)
    {
        for(k=0; k<temp2; k++)
        {
            query_temp = query_temp_buffer;

            query_temp_buffer = query[q_off +k + 1];

            sequence_temp = subject_buffer1_ptr[(s_off&0x3) + k];

            score += matrix[ (query_temp<<5) + sequence_temp ];

            if (score > maxscore)
            {
                maxscore = score;
                best_i = k;
            }

            if (score <= 0 || (maxscore - score) >= dropoff)
            {
                flag = 1;
                break;
            }

        }

        if ( (0==flag) && (n>temp2) )
        {
            for (int i = 0; i < n2; i++)
            {

                subject_buffer1=subject_buffer2;
                index=index + temp;
                subject_buffer2 = subject[index];

                sequence_temp = subject_buffer1_ptr[0];

                query_temp = query_temp_buffer;

                query_temp_buffer = query[q_off +k + 1];

                score += matrix[ query_temp*32 + sequence_temp ];

                if (score > maxscore)
                {
                    maxscore = score;
                    best_i = k;
                }

                if (score <= 0 || (maxscore - score) >= dropoff)
                {
                    flag = 1;
                    break;
                }

                k++;
                sequence_temp = subject_buffer1_ptr[1];

                query_temp = query_temp_buffer;

                query_temp_buffer = query[q_off +k + 1];

                score += matrix[ query_temp*32 + sequence_temp ];

                if (score > maxscore)
                {
                    maxscore = score;
                    best_i = k;
                }

                if (score <= 0 || (maxscore - score) >= dropoff)
                {
                    flag = 1;
                    break;
                }

                k++;
                sequence_temp = subject_buffer1_ptr[2];

                query_temp = query_temp_buffer;

                query_temp_buffer = query[q_off +k + 1];

                score += matrix[ query_temp*32 + sequence_temp ];

                if (score > maxscore)
                {
                    maxscore = score;
                    best_i = k;
                }

                if (score <= 0 || (maxscore - score) >= dropoff)
                {
                    flag = 1;
                    break;
                }

                k++;
                sequence_temp = subject_buffer1_ptr[3];

                query_temp = query_temp_buffer;

                query_temp_buffer = query[q_off +k + 1];

                score += matrix[ query_temp*32 + sequence_temp ];

                if (score > maxscore)
                {
                    maxscore = score;
                    best_i = k;
                }

                if (score <= 0 || (maxscore - score) >= dropoff)
                {
                    flag = 1;
                    break;
                }

                k++;

            }


            if(0==flag)
            {
                subject_buffer1=subject_buffer2;
                for(int j=0; j<(n1 &0x3); j++,k++)
                {
                    sequence_temp = subject_buffer1_ptr[j];

                    query_temp = query_temp_buffer;

                    query_temp_buffer = query[q_off +k + 1];

                    score += matrix[ query_temp*32 + sequence_temp ];

                    if (score > maxscore)
                    {
                        maxscore = score;
                        best_i = k;
                    }

                    if (score <= 0 || (maxscore - score) >= dropoff)
                        break;

                }
            }

        }
    }


    ExtendRightReturn ReturnValues;
    ReturnValues.length = best_i + 1;
    ReturnValues.s_last_off = s_off + k;
    ReturnValues.maxscore = maxscore;

    return ReturnValues;
}


__device__ ExtendLeftReturn s_BlastAaExtendLeft(const char* matrix,
        const int* __restrict__ subject,
        const char* __restrict__ query,
        const int s_off,
        const int q_off,
        const int dropoff,
        int maxscore)
{
    int n, best_i;
    int score = maxscore, score1, score2, score3;

    char query_temp_buffer = query[q_off], query_temp_buffer1;
    char query_temp_buffer2, query_temp_buffer3;
    int subject_buffer1, subject_buffer2;
    int temp = blockDim.x;
    int index = ((s_off)>>2)*temp;
    subject_buffer1 = subject[index];

    index -= temp;
    subject_buffer2 = subject[index ];
    char* subject_buffer1_ptr = (char*) &subject_buffer1;


    n = MIN(s_off, q_off);

    int temp2 = ((s_off)&0x3);

    best_i = n + 1;

    int k = n;

    int n1=n-temp2;

    query_temp_buffer1 = query[q_off - temp2 -2];
    query_temp_buffer2 = query[q_off - temp2 -3];
    query_temp_buffer3 = query[q_off - temp2 -4];

    int n2 = n1>>2;

    if(n>0)
    {

        char query_temp = 0,  query_temp1, flag = 0;

        int sequence_temp = 0;
        int sequence_temp1 = 0;
        int sequence_temp2 = 0;
        int sequence_temp3 = 0;

        for(int i=temp2; i>=0; i--)
        {
            query_temp = query_temp_buffer;
            query_temp_buffer = query[q_off - n+k -1];
            score += matrix[ (query_temp<<5) + subject_buffer1_ptr[i] ];

            if (score > maxscore)
            {
                maxscore = score;
                best_i = k;
            }

            if ((maxscore - score) >= dropoff)
            {
                flag = 1;
                break;
            }
            k--;
        }


        if ( (flag == 0) && (n1>0) )
        {
            for (int i = 0; i < n2; i++)
            {
                subject_buffer1=subject_buffer2;
                index=index- temp;
                subject_buffer2 = subject[index];

                sequence_temp1 = subject_buffer1&0xff0000;
                sequence_temp = subject_buffer1>>24;

                query_temp = query_temp_buffer;

                query_temp_buffer = query[q_off - n+k -4];

                sequence_temp1 >>=16;

                query_temp1 = query_temp_buffer1;

                query_temp_buffer1 = query[q_off - n+k -5];

                score += matrix[ (query_temp<<5) + sequence_temp ];
                score1 = matrix[ (query_temp1<<5) + sequence_temp1 ];

                sequence_temp2 = subject_buffer1&0xff00;
                sequence_temp3 = subject_buffer1&0x1f;

                query_temp = query_temp_buffer2;

                query_temp_buffer2 = query[q_off - n+k -6];

                score1 += score;
                sequence_temp = maxscore - score;

                sequence_temp2 >>=8;

                query_temp1 = query_temp_buffer3;

                query_temp_buffer3 = query[q_off - n+k -7];

                score2 = matrix[ (query_temp<<5) + sequence_temp2 ] ;
                score3 = matrix[ (query_temp1<<5) + sequence_temp3 ];

                if (score > maxscore)
                {
                    maxscore = score;
                    best_i = k;
                }

                if (sequence_temp >= dropoff)
                {
                    flag = 1;
                    break;
                }

                sequence_temp1 = maxscore - score1;
                score3 += score2;
                score2 += score1;

                k--;
                if (score1 > maxscore)
                {
                    maxscore = score1;
                    best_i = k;
                }

                if (sequence_temp1 >= dropoff)
                {
                    flag = 1;
                    break;
                }

                sequence_temp2 = maxscore - score2;
                score3 += score1;

                k--;
                if (score2 > maxscore)
                {
                    maxscore = score2;
                    best_i = k;
                }

                if (sequence_temp2 >= dropoff)
                {
                    flag = 1;
                    break;
                }

                sequence_temp3 = maxscore - score3;
                score = score3;

                k--;
                if (score3 > maxscore)
                {
                    maxscore = score3;
                    best_i = k;
                }

                if (sequence_temp3 >= dropoff)
                {
                    flag = 1;
                    break;
                }
                k--;

            }

            if(0==flag)
            {
                subject_buffer1 = subject_buffer2;
                for(int j=0; j<(n1 &0x3); j++)
                {
                    sequence_temp = subject_buffer1_ptr[3-j];

                    query_temp = query_temp_buffer;

                    query_temp_buffer = query[q_off - n+k -1];

                    score += matrix[ (query_temp<<5) + sequence_temp ];

                    if (score > maxscore)
                    {
                        maxscore = score;
                        best_i = k;
                    }

                    if ((maxscore - score) >= dropoff)
                        break;

                    k--;

                }
            }
        }
    }


    ExtendLeftReturn ReturnValues;
    ReturnValues.maxscore = maxscore;
    ReturnValues.length = n - best_i + 1;
    return ReturnValues;
}


