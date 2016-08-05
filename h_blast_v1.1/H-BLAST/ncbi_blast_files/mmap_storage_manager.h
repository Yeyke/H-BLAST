// includes for mmap
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include <stdlib.h>
#include <stdio.h>

#define MAX_NUM_OF_ARRAY_ELEMENTS 100

typedef struct
{
    char* start_for_mem_file;
    unsigned long size_of_mem_file_in_bytes;
    int file_id_for_mem_file;
    char mem_file_is_allocated;
} mmap_storage_manager;

static mmap_storage_manager mmap_array[MAX_NUM_OF_ARRAY_ELEMENTS];
static unsigned int no_of_mmap_used;

//static mmap_storage_manager init_mmap_storage_manager_by_info(char* addr, unsigned long size, int file_id)
//{
//    mmap_storage_manager ret_val;
//
//    return ret_val;
//
//}

static int is_mmap_allocated_by_id(unsigned int i)
{
    if(i>=no_of_mmap_used)
        return -1;
    return (mmap_array[i].mem_file_is_allocated)?1:0;
}

static void init_a_mmap_storage_manager(mmap_storage_manager* in_element)
{
    if(in_element==NULL)
    return;
    in_element->start_for_mem_file=NULL;
    in_element->size_of_mem_file_in_bytes=0;
    in_element->file_id_for_mem_file=-1;
    in_element->mem_file_is_allocated=0;
}


static void init_mmap_storage_manager_by_id(unsigned int i)
{
    if(i<MAX_NUM_OF_ARRAY_ELEMENTS)
        init_a_mmap_storage_manager(&(mmap_array[i]));
}

static void init_mmap_storage_manager_array()
{
    for(int i=0; i<MAX_NUM_OF_ARRAY_ELEMENTS; i++)
    {
        init_mmap_storage_manager_by_id(i);
    }
    no_of_mmap_used = 0;
    return;
}

static mmap_storage_manager get_mmap_storage_manager_by_id(unsigned int i)
{
    mmap_storage_manager ret;
    init_a_mmap_storage_manager(&ret);
    if(is_mmap_allocated_by_id(i)>0)
        return mmap_array[i];
    return ret;
}

static int set_mmap_storage_manager_by_id(mmap_storage_manager element, unsigned int i)
{
    if(i>=MAX_NUM_OF_ARRAY_ELEMENTS)
        return -1;
    else if(mmap_array[i].mem_file_is_allocated > 0)
        return -2;

    mmap_array[i].mem_file_is_allocated = (NULL != element.start_for_mem_file)? 1: 0;
    mmap_array[i].start_for_mem_file=element.start_for_mem_file;
    mmap_array[i].size_of_mem_file_in_bytes=element.size_of_mem_file_in_bytes;
    mmap_array[i].file_id_for_mem_file=element.file_id_for_mem_file;
    no_of_mmap_used = (no_of_mmap_used<(i+1))?(i+1):no_of_mmap_used;
    return 0;
}



static int push_a_mmap_storage_manager_to_array(mmap_storage_manager element)
{
    return set_mmap_storage_manager_by_id(element,no_of_mmap_used);
}

static int push_a_mmap_storage_manager_to_array_with_info(char* addr, unsigned long size, int file_id)
{
    mmap_storage_manager m_tmp;
    m_tmp.mem_file_is_allocated = (NULL != addr)? 1: 0;;
    m_tmp.start_for_mem_file=addr;
    m_tmp.size_of_mem_file_in_bytes=size;
    m_tmp.file_id_for_mem_file=file_id;
    return push_a_mmap_storage_manager_to_array(m_tmp);

}

static void release_mmap_storage_manager_by_id(unsigned int i)
{

    if(is_mmap_allocated_by_id(i)<1)
        return;

    int ret = -1;
    ret = munmap(mmap_array[i].start_for_mem_file, mmap_array[i].size_of_mem_file_in_bytes);
    if(-1 == ret)
    {
        fputs("munmap mmap file failed: \n", stderr);
        exit (1);
    }
    close(mmap_array[i].file_id_for_mem_file);

    init_mmap_storage_manager_by_id(i);

    return;

}

static void release_mmap_storage_manager_array()
{
    for(int i=0; i<no_of_mmap_used; i++)
        release_mmap_storage_manager_by_id(i);
    no_of_mmap_used=0;
    return;
}


