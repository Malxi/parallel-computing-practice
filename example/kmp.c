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


/*��ģʽ���������ڷ�������������Ӧ��new��newvalֵ*/
void Next(char *W,int patlen,int *nextval,pntype *pped)
{
	int i,j,plen;
   	int *next;
    
    /* �����СΪpatlen + 1�������ڴ��next */
   	if((next=(int *)malloc((patlen+1)*sizeof(int)))==NULL){
     	printf("no enough memory\n");
       	exit(1);
   	}

  	/*����next��nextval*/
    /* ��ʼ��next[0]��nextval[0]Ϊ -1 */
    next[0]=nextval[0]=-1;
  	j=1;
    /* ����W��ÿ��λ��j������next���� */
   	while(j<=patlen){
   		/* i���ȸ�ֵΪnext[j-1], ǰһ������next */
        i=next[j-1];
        /* Ѱ��iʹ��W[i] �� w[j-1] ��� ���� Ϊ -1 */
      	while(i!=(-1) && W[i]!=W[j-1]) i=next[i];
      	/* ����j��nextΪi+1 */
        next[j]=i+1;

      	if(j!=patlen){
            /* ��jλ���ַ���i+1λ�ò���ȣ���ֵΪi+1 */
          	if( W[j]!=W[i+1])
              	nextval[j]=i+1;
           	else
            /* ��jλ���ַ���i+1λ����ȣ���ֵΪnextval[i+1] */
    			nextval[j]=nextval[i+1]; 
      	}
        /* ��һ��λ�� */
        j++;
  	} 
    
    /* ����W������ַ�������д��pped */
    pped->pedlen=patlen-next[patlen]; 
 	pped->pednum=(int)(patlen/pped->pedlen); 
   	pped->psuffixlen=patlen%pped->pedlen; 

  	free(next);
}

/*�Ľ���KMP�㷨*/
void kmp(char *T,char*W,int textlen,int patlen,int *nextval,pntype *pped,int prefix_flag,int matched_num,int *match,int *prefixlen)
{
	int i,j;
    
    /* ƥ����Ŀ */
  	i=matched_num;              
  	j=matched_num;           
    
    /* iС��textlenʱ����ƥ�� */
    while(i<textlen)
    {
    	/* ǰ׺��־Ϊ1������ ������-1 �����ı�ʣ������������ƥ�䣬��Ҫ��һ������������ */
        if((prefix_flag==1)&&((patlen-j)>(textlen-i))) {
          	break;
		}
        
        /* Ѱ��j ʹ���ı�iλ����ƥ�䴮��jλ����� */
        while(j!=(-1) && W[j]!=T[i])  j=nextval[j];
        
        /* j = patlen - 1 ʱƥ��ɹ� */
        if(j==(patlen-1)) {  
            /* ��match�����ƥ��ɹ�����ʼλ�ñ��Ϊ1 */
			match[i-(patlen-1)]=1;
            /* ƥ�䴮����������ƥ�䴮��ǰ׺��Ϊ1�����j����Ϊ-1����������Ϊ���һ��������ʼpatlen-1-pped->pedlen */
        	if(pped->pednum+pped->psuffixlen==1)
                j = -1; 
            else                                   
                j=patlen-1-pped->pedlen; 
        }
   		j++;
      	i++; 
   	}
    /* ����ǰ׺����Ϊ j ,������ƥ�����j = 0������ǰ׺�ı�׼�£�����ǰ����ʣ����ַ���Ŀ����ƥ����ɣ���������һ������������ */
   	(*prefixlen)=j;
}

/*�ع�ģʽ���Լ�next����*/
void Rebuild_info(int patlen,pntype *pped,int *nextval,char *W) 
{ 
	int i; 
   	if (pped->pednum == 1)
        /* ֻ��һ�����ڣ���ĩǰ׺��������һ����Χĩβ���ع�ƥ�䴮�ɹ� */
   		memcpy(W+pped->pedlen,W,pped->psuffixlen); 
	else {  
       	/* ��ʱ�������������ڣ��ع��ڶ������� */
        memcpy(W+pped->pedlen,W,pped->pedlen);
       	
        /* ����ʣ�µ����ڣ���Ҫ����nextval���飬��Ϊ֮ǰ��nextvalΪ�������ڳ��� */
        for (i=3; i<=pped->pednum; i++){ 
        	memcpy(W+(i-1)*pped->pedlen,W,pped->pedlen);
           	memcpy(nextval+(i-1)*pped->pedlen,nextval+pped->pedlen,pped->pedlen*sizeof(int));
      	} 
        /* �ع�ʣ�µ�ĩβǰ׺ */
       	if(pped->psuffixlen!=0){
       		memcpy(W+(i-1)*pped->pedlen,W,pped->psuffixlen);
           	memcpy(nextval+(i-1)*pped->pedlen,nextval+pped->pedlen,pped->psuffixlen*sizeof(int));            
       	}
 	} 
} 
 
/*�����ı���*/
void gen_string(int strlen,int pedlen,char *string,int seed)
{
	int suffixlen,num,i,j;

   	srand(seed);
    /* ���ɳ���Ϊpedlen���ַ������浽��ʼλ�� */
   	for(i=0;i<pedlen;i++){
    	num=rand()%26;
        string[i]='a'+num;
   	}
    /* ����ʼ��pedlen�����ַ�����������ÿ��pedlen����λ�� */
   	for(j=1;j<(int)(strlen/pedlen);j++)
    	strncpy(string+j*pedlen,string,pedlen);
   	/* ��pedlen���ȵ��ַ���ǰ׺��ȫ������ַ�λ�� */
    if((suffixlen=strlen%pedlen)!=0)
    	strncpy(string+j*pedlen,string,suffixlen);
}  

/*���ļ�����ģʽ����Ϣ*/ 
void GetFile(char *filename,char **place, int *number) 
{ 
	FILE *fp;
    struct stat statbuf;  

    if ((fp=fopen(filename,"rb"))==NULL) {
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 
    
    /* ��ȡ�ļ���Ϣ */
    fstat(fileno(fp),&statbuf);
    
    /* �����ļ���С���ڴ� */
    if(((*place)=(char *)malloc(sizeof(char)*statbuf.st_size)) == NULL){
		printf ("Error alloc memory\n");
        exit(1);
 	}
    
    /* ��ȡ�ļ����ŵ�place */
   	if (fread((*place),1,statbuf.st_size,fp)!=statbuf.st_size){
		printf ("Error in reading num\n"); 
        exit(0); 
	} 
    fclose (fp); 
    /* �޸�number��ֵΪ�ļ����� */
    (*number)=statbuf.st_size; 
} 

/*��ӡ���в�����Ϣ*/
void PrintFile_info(char *filename,char *T,int id)
{ 
	FILE *fp; 
	int i;
    
	if ((fp=fopen(filename,"a"))==NULL){
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 
    
    /* ��ÿ��Node�ϵ��ַ���ӡ���� */
 	fprintf (fp,"The Text on node %d is %s .\n",id,T); 
   	
	fclose (fp); 
} 


/*��ӡƥ����*/
void PrintFile_res(char *filename,int *t,int len,int init,int id)
{ 
	FILE *fp; 
	int i;

	if ((fp=fopen(filename,"a"))==NULL){
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 
    /* ��ƥ��ɹ���λ��������ļ� */
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
   	/* ��ȡ�ı������� */
   	textlen=atoi(argv[1]);
    /* ���ı����Ⱦ��ָ�ÿ������ */
   	textlen=textlen/numprocs;
    /* ��ȡ�ı������ڳ��� */
  	pedlen=atoi(argv[2]);
    /* �����СΪtextlen���ڴ� */
   	if((T=(char *)malloc(textlen*sizeof(char)))==NULL){
     	printf("no enough memory\n");
       	exit(1);
   	}
    
    /* ���ɳ���Ϊtextlen���ı�������Ϊpedlen������ʹ�õ�ǰ���̵�ID�������ȷ���ģ��������T�� */
   	gen_string(textlen,pedlen,T,myid);
    
    /* ʵ������0~(n-1)��ӡ���ɵ��ַ��ı���Ϣ */
	if(myid==0){
        /* ����0�Ž��̴�ӡ���ɵ��ַ��ı���Ϣ */
		PrintFile_info("match_result",T,myid);
		/* ������������1����ready���͸�1�Ž��̣�tag����Ϊ0 */
        if(numprocs>1)
			MPI_Send(&ready,1,MPI_INT,1,0,MPI_COMM_WORLD);
	}
   	else{
        /* ��0�Ž��̽��������Ա����ǰһ�����̵�ready��Ϣ */
  		MPI_Recv(&ready,1,MPI_INT,myid-1,myid-1,MPI_COMM_WORLD,&status);
        /* ���̴�ӡ���ɵ��ַ��ı���Ϣ */
		PrintFile_info("match_result",T,myid);
        /* ���������һ��������ready��Ϣ���͸���һ������myid+1��tagΪmyid */
		if(myid!=numprocs-1)
			MPI_Send(&ready,1,MPI_INT,myid+1,myid,MPI_COMM_WORLD);
	}

	printf("\n");
   	
    /* �����СΪtextlen�������ڴ棬����ֵ��match */
	if((match=(int *)malloc(textlen*sizeof(int)))==NULL){
		printf("no enough memory\n");
		exit(1);
	}
  
  	/*������0����ģʽ��������¼���в���*/
   	if(myid==0){
        /* ��ӡ���������ܵ��ı������� */
  		printf("processor num = %d \n",numprocs);
    	printf("textlen = %d\n",textlen*numprocs); 
        
        /* ��ȡpattern.bat���ļ����������ݷŵ�W���ļ�����Ϊpatlen */ 
        GetFile("pattern.dat",&W,&patlen); 
    	/* ��ӡpatlen�ĳ��� */
        printf("patlen= %d\n",patlen); 
        
        /* �����СΪpatlen���ڴ��nextval�����ڴ洢nextֵ */
    	if((nextval=(int *)malloc(patlen*sizeof(int)))==NULL){
        	printf("no enough memory\n");
           	exit(1);
        }
		/*��ģʽ�����з�������Ӧ���㷨14.6���裨1��*/
      	Next(W,patlen,nextval,&pped);
        if(numprocs>1){
            /* ���봮������Ϊ1 */
        	if (pped.pednum==1) 
           		nextlen_send = patlen;
            else 
            /* ���봮������Ϊ1, nextlen_send�Ĵ�СΪ���봮���ڵ����� */
        		nextlen_send = pped.pedlen*2;
        }
    }

	/*���������������ģʽ������Ϣ����Ӧ���㷨14.6���裨2��*/
  	if(numprocs>1){
        /* 0�Ž��������н��̹㲥����ƥ�䴮���� *//
     	MPI_Bcast(&patlen, 1, MPI_INT, 0, MPI_COMM_WORLD);  
        
        /* ����0�Ž��̷ֱ������СΪpatlen�������ڴ��nextval��W */
  		if(myid!=0)
    		if(((nextval=(int *)malloc(patlen*sizeof(int)))==NULL)
				||((W=(char *)malloc(patlen*sizeof(char)))==NULL)){
           		printf("no enough memory\n");
            	exit(1);
            }
        
        /* ͬ�����н��� */
 	 	MPI_Barrier(MPI_COMM_WORLD);
        /* 0�Ž��̹㲥pped��Ϣ */
    	MPI_Bcast(&pped,3,MPI_INT,0,MPI_COMM_WORLD);
        /* 0�Ž��̹㲥nextlen��Ϣ */
    	MPI_Bcast(&nextlen_send,1,MPI_INT,0,MPI_COMM_WORLD);
        /* 0�Ž��̹㲥nextval��Ϣ */
    	MPI_Bcast(nextval,nextlen_send,MPI_INT,0,MPI_COMM_WORLD); 
        /* 0�Ž��̹㲥W��Ϣ, ֻ�����봮�ĵ�һ�����ڲ��� */
    	MPI_Bcast(W,pped.pedlen,MPI_CHAR,0,MPI_COMM_WORLD);
   	}
    
    /* ͬ�����н��� */
    MPI_Barrier(MPI_COMM_WORLD);

   	/*�����޸Ĺ���KMP�㷨���оֲ���ƥ�䣬��Ӧ���㷨14.6���裨3��*/
  	if(numprocs==1) {
        /* ����һ�����̣� ����KMP */
  		kmp(T,W,textlen,patlen,nextval,&pped,1,0,match+patlen-1,&prefixlen);
   	}
   	else { 
    	if(myid!=0)
    		/*�����������ֱ���ݲ��ִ������Լ�������Ϣ�ع�ģʽ��*/
        	Rebuild_info(patlen,&pped,nextval,W); 
    	if(myid!=numprocs-1)
            /* �������һ�����̵�ǰ׺flag = 0 */
  			kmp(T,W,textlen,patlen,nextval,&pped,0,0,match+patlen-1,&prefixlen);
		else
  			kmp(T,W,textlen,patlen,nextval,&pped,1,0,match+patlen-1,&prefixlen);

   		/* ͬ�����н��� */
        MPI_Barrier(MPI_COMM_WORLD);

		/*�������������жμ�ƥ�䣬��Ӧ���㷨14.6���裨4��*/
    	if(myid<numprocs-1) 
            /* ��ƥ���ǰ׺���͵���һ�����̽��д��� */
        	MPI_Send(&prefixlen,1,MPI_INT,myid+1,99,MPI_COMM_WORLD); 

    	if(myid>0)
            /* ��Ŵ���0�Ľ��̽���֮ǰ���̵�ǰ׺ */
    		MPI_Recv(&prefixlen,1,MPI_INT,myid-1,99,MPI_COMM_WORLD,&status); 

    	/* ͬ�����н��� */
        MPI_Barrier(MPI_COMM_WORLD);

    	/* �����̱�Ŵ���0������֮ǰ���̷��͵�ǰ׺��������0������Ҫ����ƥ�䣬��ʱ�Ѿ�ƥ�����ĿΪprefixlen */
        if((myid>0)&&(prefixlen!=0))  
   			kmp(T-prefixlen,W,prefixlen+patlen-1,patlen,nextval,&pped,1,prefixlen,match+patlen-1-prefixlen,&prefixlen); 
        
        /* ͬ�����н��� */
   		MPI_Barrier(MPI_COMM_WORLD);
   	}

	/*���ƥ����*/
	if(myid==0){
		PrintFile_res("match_result",match+patlen-1,textlen-patlen+1,0,myid);
		if(numprocs>1)
            /* ֪ͨ1�Ž��̴�ӡ��� */
			MPI_Send(&ready,1,MPI_INT,1,0,MPI_COMM_WORLD);
	}
   	else{
        /* ��0������Ҫ�ȴ�ǰ����̵��źţ�Ȼ���ӡ */
  		MPI_Recv(&ready,1,MPI_INT,myid-1,myid-1,MPI_COMM_WORLD,&status);
		PrintFile_res("match_result",match,textlen,myid*textlen-patlen+1,myid);
		if(myid!=numprocs-1)
            /* �������Ľ�����Ҫ֪ͨ��һ������ */
			MPI_Send(&ready,1,MPI_INT,myid+1,myid,MPI_COMM_WORLD);
	}

	free(T);
    free(W);
    free(nextval);
    MPI_Finalize(); 
 } 
