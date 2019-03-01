47554.console  
#include <malloc.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#define  MAX(m,n)    (m>n?m:n)

typedef struct{
  	int pedlen;
  	int psuffixlen;
  	int pednum;
}pntype;


/*对模式串进行周期分析，并计算相应的new和newval值*/
void Next(char *W,int patlen,int *nextval,pntype *pped)
{
	int i,j,plen;
   	int *next;
    
    /* 分配大小为patlen + 1的整型内存给next */
   	if((next=(int *)malloc((patlen+1)*sizeof(int)))==NULL){
     	printf("no enough memory\n");
       	exit(1);
   	}

  	/*计算next和nextval*/
    /* 初始化next[0]和nextval[0]为 -1 */
    next[0]=nextval[0]=-1;
  	j=1;
    /* 对于W的每个位置j都进行next计算 */
   	while(j<=patlen){
   		/* i首先赋值为next[j-1], 前一个数的next */
        i=next[j-1];
        /* 寻找i使得W[i] 与 w[j-1] 相等 或者 为 -1 */
      	while(i!=(-1) && W[i]!=W[j-1]) i=next[i];
      	/* 设置j的next为i+1 */
        next[j]=i+1;

      	if(j!=patlen){
            /* 若j位置字符和i+1位置不想等，则赋值为i+1 */
          	if( W[j]!=W[i+1])
              	nextval[j]=i+1;
           	else
            /* 若j位置字符和i+1位置相等，则赋值为nextval[i+1] */
    			nextval[j]=nextval[i+1]; 
      	}
        /* 下一个位置 */
        j++;
  	} 
    
    /* 将对W输入的字符串分析写入pped */
    pped->pedlen=patlen-next[patlen]; 
 	pped->pednum=(int)(patlen/pped->pedlen); 
   	pped->psuffixlen=patlen%pped->pedlen; 

  	free(next);
}

/*改进的KMP算法*/
void kmp(char *T,char*W,int textlen,int patlen,int *nextval,pntype *pped,int prefix_flag,int matched_num,int *match,int *prefixlen)
{
	int i,j;
    
    /* 匹配数目 */
  	i=matched_num;              
  	j=matched_num;           
    
    /* i小于textlen时进行匹配 */
    while(i<textlen)
    {
    	/* 前缀标志为1，并且 周期数-1 大于文本剩余数量，结束匹配，需要下一个进程来处理 */
        if((prefix_flag==1)&&((patlen-j)>(textlen-i))) {
          	break;
		}
        
        /* 寻找j 使得文本i位置与匹配串的j位置相等 */
        while(j!=(-1) && W[j]!=T[i])  j=nextval[j];
        
        /* j = patlen - 1 时匹配成功 */
        if(j==(patlen-1)) {  
            /* 将match数组的匹配成功串起始位置标记为1 */
			match[i-(patlen-1)]=1;
            /* 匹配串的周期数和匹配串的前缀和为1，则把j设置为-1，否则设置为最后一个周期起始patlen-1-pped->pedlen */
        	if(pped->pednum+pped->psuffixlen==1)
                j = -1; 
            else                                   
                j=patlen-1-pped->pedlen; 
        }
   		j++;
      	i++; 
   	}
    /* 设置前缀长度为 j ,正常地匹配完成j = 0，在有前缀的标准下，若当前进程剩余的字符数目不足匹配完成，则需由下一个进程来处理 */
   	(*prefixlen)=j;
}

/*重构模式串以及next函数*/
void Rebuild_info(int patlen,pntype *pped,int *nextval,char *W) 
{ 
	int i; 
   	if (pped->pednum == 1)
        /* 只有一个周期，则将末前缀拷贝至第一个周围末尾，重构匹配串成功 */
   		memcpy(W+pped->pedlen,W,pped->psuffixlen); 
	else {  
       	/* 此时至少有两个周期，重构第二个周期 */
        memcpy(W+pped->pedlen,W,pped->pedlen);
       	
        /* 对于剩下的周期，需要拷贝nextval数组，因为之前的nextval为两倍周期长度 */
        for (i=3; i<=pped->pednum; i++){ 
        	memcpy(W+(i-1)*pped->pedlen,W,pped->pedlen);
           	memcpy(nextval+(i-1)*pped->pedlen,nextval+pped->pedlen,pped->pedlen*sizeof(int));
      	} 
        /* 重构剩下的末尾前缀 */
       	if(pped->psuffixlen!=0){
       		memcpy(W+(i-1)*pped->pedlen,W,pped->psuffixlen);
           	memcpy(nextval+(i-1)*pped->pedlen,nextval+pped->pedlen,pped->psuffixlen*sizeof(int));            
       	}
 	} 
} 
 
/*生成文本串*/
void gen_string(int strlen,int pedlen,char *string,int seed)
{
	int suffixlen,num,i,j;

   	srand(seed);
    /* 生成长度为pedlen的字符串，存到起始位置 */
   	for(i=0;i<pedlen;i++){
    	num=rand()%26;
        string[i]='a'+num;
   	}
    /* 将起始的pedlen长度字符串，拷贝到每个pedlen周期位置 */
   	for(j=1;j<(int)(strlen/pedlen);j++)
    	strncpy(string+j*pedlen,string,pedlen);
   	/* 用pedlen长度的字符串前缀补全多余的字符位置 */
    if((suffixlen=strlen%pedlen)!=0)
    	strncpy(string+j*pedlen,string,suffixlen);
}  

/*从文件读入模式串信息*/ 
void GetFile(char *filename,char **place, int *number) 
{ 
	FILE *fp;
    struct stat statbuf;  

    if ((fp=fopen(filename,"rb"))==NULL) {
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 
    
    /* 读取文件信息 */
    fstat(fileno(fp),&statbuf);
    
    /* 分配文件大小的内存 */
    if(((*place)=(char *)malloc(sizeof(char)*statbuf.st_size)) == NULL){
		printf ("Error alloc memory\n");
        exit(1);
 	}
    
    /* 读取文件，放到place */
   	if (fread((*place),1,statbuf.st_size,fp)!=statbuf.st_size){
		printf ("Error in reading num\n"); 
        exit(0); 
	} 
    fclose (fp); 
    /* 修改number的值为文件长度 */
    (*number)=statbuf.st_size; 
} 

/*打印运行参数信息*/
void PrintFile_info(char *filename,char *T,int id)
{ 
	FILE *fp; 
	int i;
    
	if ((fp=fopen(filename,"a"))==NULL){
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 
    
    /* 将每个Node上的字符打印出来 */
 	fprintf (fp,"The Text on node %d is %s .\n",id,T); 
   	
	fclose (fp); 
} 


/*打印匹配结果*/
void PrintFile_res(char *filename,int *t,int len,int init,int id)
{ 
	FILE *fp; 
	int i;

	if ((fp=fopen(filename,"a"))==NULL){
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 
    /* 将匹配成功的位置输出到文件 */
 	fprintf (fp,"This is the match result on node %d \n",id); 
   	for (i=0; i<=len-1; i++) 
    	if(t[i]==1)
 			fprintf (fp,"(%d)  +\n",i+init); 
    	else
  			fprintf (fp,"(%d)  -\n",i+init);
	fclose (fp); 
} 

void main(int argc,char *argv[]) 
{ 
	char *T,*W; 
	int	*nextval,*match; 
  	int	textlen,patlen,pedlen,nextlen_send; 
   	pntype pped; 
 	int	i,myid,numprocs,prefixlen,ready; 
  	MPI_Status  status; 

	MPI_Init(&argc,&argv); 
   	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
   	MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

   	nextlen_send=0;
   	ready=1;
   	/* 读取文本串长度 */
   	textlen=atoi(argv[1]);
    /* 将文本长度均分给每个进程 */
   	textlen=textlen/numprocs;
    /* 读取文本串周期长度 */
  	pedlen=atoi(argv[2]);
    /* 分配大小为textlen的内存 */
   	if((T=(char *)malloc(textlen*sizeof(char)))==NULL){
     	printf("no enough memory\n");
       	exit(1);
   	}
    
    /* 生成长度为textlen，文本串周期为pedlen，种子使用当前进程的ID，结果是确定的，结果存入T中 */
   	gen_string(textlen,pedlen,T,myid);
    
    /* 实现有序将0~(n-1)打印生成的字符文本信息 */
	if(myid==0){
        /* 对于0号进程打印生成的字符文本信息 */
		PrintFile_info("match_result",T,myid);
		/* 若进程数多余1，则将ready发送给1号进程，tag设置为0 */
        if(numprocs>1)
			MPI_Send(&ready,1,MPI_INT,1,0,MPI_COMM_WORLD);
	}
   	else{
        /* 非0号进程接收来自自编号是前一个进程的ready消息 */
  		MPI_Recv(&ready,1,MPI_INT,myid-1,myid-1,MPI_COMM_WORLD,&status);
        /* 进程打印生成的字符文本信息 */
		PrintFile_info("match_result",T,myid);
        /* 若不是最后一个进程则将ready信息发送给下一个进程myid+1，tag为myid */
		if(myid!=numprocs-1)
			MPI_Send(&ready,1,MPI_INT,myid+1,myid,MPI_COMM_WORLD);
	}

	printf("\n");
   	
    /* 申请大小为textlen的整型内存，并赋值给match */
	if((match=(int *)malloc(textlen*sizeof(int)))==NULL){
		printf("no enough memory\n");
		exit(1);
	}
  
  	/*处理器0读入模式串，并记录运行参数*/
   	if(myid==0){
        /* 打印进程数及总的文本串长度 */
  		printf("processor num = %d \n",numprocs);
    	printf("textlen = %d\n",textlen*numprocs); 
        
        /* 读取pattern.bat的文件，将其内容放到W，文件长度为patlen */ 
        GetFile("pattern.dat",&W,&patlen); 
    	/* 打印patlen的长度 */
        printf("patlen= %d\n",patlen); 
        
        /* 分配大小为patlen的内存给nextval，用于存储next值 */
    	if((nextval=(int *)malloc(patlen*sizeof(int)))==NULL){
        	printf("no enough memory\n");
           	exit(1);
        }
		/*对模式串进行分析，对应于算法14.6步骤（1）*/
      	Next(W,patlen,nextval,&pped);
        if(numprocs>1){
            /* 输入串的周期为1 */
        	if (pped.pednum==1) 
           		nextlen_send = patlen;
            else 
            /* 输入串的周期为1, nextlen_send的大小为输入串周期的两倍 */
        		nextlen_send = pped.pedlen*2;
        }
    }

	/*向各个处理器播送模式串的信息，对应于算法14.6步骤（2）*/
  	if(numprocs>1){
        /* 0号进程向所有进程广播输入匹配串长度 *//
     	MPI_Bcast(&patlen, 1, MPI_INT, 0, MPI_COMM_WORLD);  
        
        /* 非零0号进程分别申请大小为patlen的整型内存给nextval和W */
  		if(myid!=0)
    		if(((nextval=(int *)malloc(patlen*sizeof(int)))==NULL)
				||((W=(char *)malloc(patlen*sizeof(char)))==NULL)){
           		printf("no enough memory\n");
            	exit(1);
            }
        
        /* 同步所有进程 */
 	 	MPI_Barrier(MPI_COMM_WORLD);
        /* 0号进程广播pped信息 */
    	MPI_Bcast(&pped,3,MPI_INT,0,MPI_COMM_WORLD);
        /* 0号进程广播nextlen信息 */
    	MPI_Bcast(&nextlen_send,1,MPI_INT,0,MPI_COMM_WORLD);
        /* 0号进程广播nextval信息 */
    	MPI_Bcast(nextval,nextlen_send,MPI_INT,0,MPI_COMM_WORLD); 
        /* 0号进程广播W信息, 只有输入串的第一个周期部分 */
    	MPI_Bcast(W,pped.pedlen,MPI_CHAR,0,MPI_COMM_WORLD);
   	}
    
    /* 同步所有进程 */
    MPI_Barrier(MPI_COMM_WORLD);

   	/*调用修改过的KMP算法进行局部串匹配，对应于算法14.6步骤（3）*/
  	if(numprocs==1) {
        /* 对于一个进程， 调用KMP */
  		kmp(T,W,textlen,patlen,nextval,&pped,1,0,match+patlen-1,&prefixlen);
   	}
   	else { 
    	if(myid!=0)
    		/*各个处理器分别根据部分串数据以及周期信息重构模式串*/
        	Rebuild_info(patlen,&pped,nextval,W); 
    	if(myid!=numprocs-1)
            /* 不是最后一个进程的前缀flag = 0 */
  			kmp(T,W,textlen,patlen,nextval,&pped,0,0,match+patlen-1,&prefixlen);
		else
  			kmp(T,W,textlen,patlen,nextval,&pped,1,0,match+patlen-1,&prefixlen);

   		/* 同步所有进程 */
        MPI_Barrier(MPI_COMM_WORLD);

		/*各个处理器进行段间匹配，对应于算法14.6步骤（4）*/
    	if(myid<numprocs-1) 
            /* 将匹配的前缀发送到下一个进程进行处理 */
        	MPI_Send(&prefixlen,1,MPI_INT,myid+1,99,MPI_COMM_WORLD); 

    	if(myid>0)
            /* 编号大于0的进程接收之前进程的前缀 */
    		MPI_Recv(&prefixlen,1,MPI_INT,myid-1,99,MPI_COMM_WORLD,&status); 

    	/* 同步所有进程 */
        MPI_Barrier(MPI_COMM_WORLD);

    	/* 若进程编号大于0，并且之前进程发送的前缀数不等于0，则需要进行匹配，此时已经匹配的数目为prefixlen */
        if((myid>0)&&(prefixlen!=0))  
   			kmp(T-prefixlen,W,prefixlen+patlen-1,patlen,nextval,&pped,1,prefixlen,match+patlen-1-prefixlen,&prefixlen); 
        
        /* 同步所有进程 */
   		MPI_Barrier(MPI_COMM_WORLD);
   	}

	/*输出匹配结果*/
	if(myid==0){
		PrintFile_res("match_result",match+patlen-1,textlen-patlen+1,0,myid);
		if(numprocs>1)
            /* 通知1号进程打印完毕 */
			MPI_Send(&ready,1,MPI_INT,1,0,MPI_COMM_WORLD);
	}
   	else{
        /* 非0进程需要等待前面进程的信号，然后打印 */
  		MPI_Recv(&ready,1,MPI_INT,myid-1,myid-1,MPI_COMM_WORLD,&status);
		PrintFile_res("match_result",match,textlen,myid*textlen-patlen+1,myid);
		if(myid!=numprocs-1)
            /* 不是最后的进程需要通知下一个进程 */
			MPI_Send(&ready,1,MPI_INT,myid+1,myid,MPI_COMM_WORLD);
	}

	free(T);
    free(W);
    free(nextval);
    MPI_Finalize(); 
 } 
