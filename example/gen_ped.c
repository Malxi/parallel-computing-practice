#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

int main(int argc,char *argv[])
{
	int strlen,pedlen,suffixlen,num,i,j;
   	char *string;
   	FILE *fp;
    
    /* 获取文本串长度，文本串周期长度，随机数种子 */
   	strlen=atoi(argv[1]);
   	pedlen=atoi(argv[2]);
   	srand(atoi(argv[3]));
    
    /* 申请大小strlen的字符内存 */
   	string=(char*)malloc(strlen*sizeof(char));
   	if(string==NULL){
      	printf("malloc error\n");
      	exit(1);
   	}
    
    /* 在文本起始生成长度为pedlen的字符串 */
   	for(i=0;i<pedlen;i++){
        num=rand()%26;
        string[i]='a'+num;
  	}
    
    /* 将起始的字符串拷贝到每个pedlen周期的位置 */
  	for(j=1;j<(int)(strlen/pedlen);j++)
       	strncpy(string+j*pedlen,string,pedlen);
    
    /* 若文本长度模周期长度不为0，在最后一个周期之后拷贝部分pedlen字符串的前缀 */
   	if((suffixlen=strlen%pedlen)!=0)
   		strncpy(string+j*pedlen,string,suffixlen);
    
    /* 写入到文件 */
  	if((fp=fopen(argv[4],"w"))!=NULL){
     	fprintf(fp,"%s",string); 
     	fclose(fp);
   	}
   	else{
     	printf("file open error\n");
     	exit(1);
   	}

	return;
}
