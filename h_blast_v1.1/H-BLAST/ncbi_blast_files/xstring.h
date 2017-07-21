// includes for string processing
#include <string.h>

char* blastp_app_path;

typedef struct{
    char* str;
    int str_leng;
}xstring;

void init_xstring(xstring* in)
{
    if(NULL==in)
        return;
    in->str_leng=0;
    in->str=(char*)calloc(1,1);
    in->str[0]='\0';

}

void free_xstring(xstring* in)
{
        if(NULL==in)
        return;
        if(NULL!=in->str)
            free(in->str);
        in->str=NULL;
        in->str_leng=0;

}

xstring get_dir(char*  in)
{
    char* i=strrchr(in, '/');
    xstring out;

    init_xstring(&out);

    if(i!=NULL)
    {

//        out.str = (char*)calloc((size_t)(i-in+1), sizeof(char));
        out.str = (char*)realloc(out.str, (size_t)(i-in+2)*sizeof(char));
        if(NULL == out.str)
        {
            printf("Error, can not allocate memory for out string\n");
            exit(-1);
        }
        out.str_leng= (int)(i-in+2);
        memcpy(out.str, in, (int)(i-in+1));
        out.str[i-in+1]='\0';
    }

     return out;
}

xstring get_filename_only(char*  in, int length)
{
    char* i=strrchr(in, '/');
    xstring out;

    init_xstring(&out);

    if(i!=NULL)
    {
        out.str = (char*)calloc((size_t)(length-(i-in+1)), sizeof(char));
        if(NULL == out.str)
        {
            printf("Error, can not allocate memory for out string\n");
            exit(-1);
        }
        out.str_leng= length-(i-in+1);
        memcpy(out.str, i+1, length-(i-in+1));
    }
    else
    {
        out.str = (char*)calloc((size_t)(length), sizeof(char));
        if(NULL == out.str)
        {
            printf("Error, can not allocate memory for out string\n");
            exit(-1);
        }
        out.str_leng= length;
        memcpy(out.str, in, length);
    }

     return out;
}

void chop_last(char* in)
{
   if(NULL==in)
	return;
   if(strlen(in)>0)
	in[strlen(in)-1]='\0';
   return;
}

char* xstrcat(const xstring in1, const char* in2)
{

    if( (NULL==in1.str) && (NULL==in2))
    {
        printf("the string is empty\n");
        exit(-1);
    }

    xstring out;
    out.str = (char*)calloc(strlen(in2)+strlen(in1.str)+1,1);
    memcpy(out.str, in1.str, strlen(in1.str));
    strcat(out.str, in2);

    return out.str;
}
