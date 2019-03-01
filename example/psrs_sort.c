47417.console

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define INIT_TYPE 10
#define ALLTOONE_TYPE 100
#define ONETOALL_TYPE 200
#define MULTI_TYPE 300
#define RESULT_TYPE 400
#define RESULT_LEN 500
#define MULTI_LEN 600

int Spt;
long DataSize;
int *arr,*arr1;
int mylength;
int *index;
int *temp1;

main(int argc,char* argv[])
{
    long BaseNum = 1;
    int PlusNum;
    int MyID;

    MPI_Init(&argc,&argv);
    /* 获取当前进程编号 */
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    
    /* PlusNum 为60,数据大小为BaseNum*PlusNum */
    PlusNum=60;
    DataSize = BaseNum*PlusNum;
    
    /* 零号进程输出 PlusNum */
    if (MyID==0)
        printf("The DataSize is : %lu\n",PlusNum);
    
    /* Psrs_Main调用 */
    Psrs_Main();
    
    /* 零号进程输出\n */
    if (MyID==0)
        printf("\n");

    MPI_Finalize();
}


Psrs_Main( )
{
    int i,j;
    int MyID,SumID;
    int n,c1,c2,c3,c4,k,l;
    FILE * fp;
    int ready;
    MPI_Status status[32*32*2];
    MPI_Request request[32*32*2];
    
    /* MyId是当前进程编号,SumID是进程数目 */
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Comm_size(MPI_COMM_WORLD,&SumID);
    
    /* Spt为最后一个进程编号 */
    Spt = SumID-1;

	/*初始化参数*/
    /* 分配两份DataSize大小内存,arr指向开始,arr1指向第二份开始 */
    arr = (int *)malloc(2*DataSize*sizeof(int));
    if (arr==0) merror("malloc memory for arr error!");
    arr1 = &arr[DataSize];
    
    /* 若进程数大于1,则分配SumID*Spt内存给temp1,2*SumID给index, temp1存储代表元素大小至少为SumID*Spt， index存段首未，大小为2*SumID */
    if (SumID>1)
    {
        temp1 = (int *)malloc(sizeof(int)*SumID*Spt);
        if (temp1==0) merror("malloc memory for temp1 error!");
        index = (int *)malloc(sizeof(int)*2*SumID);
        if (index==0) merror("malloc memory for index error!");
    }
    
    /* 等待所有进程都将变量初始化,内存分配完成 */
    MPI_Barrier( MPI_COMM_WORLD);
    
    /* 计算当前进程排序长度,设置随机数种子为当前进程编号 */
    mylength = DataSize / SumID;
    srand(MyID);
    
    /* 输出当前节点及该节点的输入数据arr,由随机数产生,由于种子是固定的,在节点数相同时,输入是相同的 */
    printf("This is node %d \n",MyID);
    printf("On node %d the input data is:\n",MyID);
    for (i=0;i<mylength;i++)
    {
        arr[i] = (int)rand();
        printf("%d : ",arr[i]);
    }
    printf("\n");

	/*每个处理器将自己的n/P个数据用串行快速排序(Quicksort)，得到一个排好序的序列，对应于算法13.5步骤（1）*/
    MPI_Barrier( MPI_COMM_WORLD);
    quicksort(arr,0,mylength - 1);
    MPI_Barrier( MPI_COMM_WORLD);

	/*每个处理器从排好序的序列中选取第w，2w，3w，…，(P-1)w个共P-1个数据作为代表元素，其中w=n/P*P，对应于算法13.5步骤（2）*/
    if (SumID>1)
    {
        /* 等待所有进程都将自己的数据排序完成 */
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* n为抽样的间隔，每个进程都抽取spt(p-1)个代表元素，并将其放到temp1 */
        n = (int)(mylength/(Spt+1));
        for (i=0;i<Spt;i++)
            temp1[i] = arr[(i+1)*n-1];

        /* 等待所有进程抽样完毕 */
        MPI_Barrier(MPI_COMM_WORLD);

        if (MyID==0)
        {
			/* 每个处理器将选好的代表元素送到处理器P0中，对应于算法13.5步骤（3） */
            /* Irecv不需要等待发送或接收消息完成就可以执行其他任务，操作放置到request请求句柄中，每个接受一个局部 */
            /* 接收数据放置在temp1的每个spt偏移位置 */
            j = 0;
            for (i=1;i<SumID;i++)
                MPI_Irecv(&temp1[i*Spt], sizeof(int)*Spt, MPI_CHAR, i,ALLTOONE_TYPE+i, MPI_COMM_WORLD, &request[j++]);
            /* 等待接收其余进程的代表元素 */
            MPI_Waitall(SumID-1, request, status);

			/* 处理器P0将上一步送来的P段有序的数据序列做P路归并，再选择排序后的第P-1，2(P-1)，…，(P-1)(P-1)个共P-1个主元，，对应于算法13.5步骤（3）*/
            /* 等待其余进程到达该位置 */
            MPI_Barrier(MPI_COMM_WORLD);
            quicksort(temp1,0,SumID*Spt-1);
            /* 同步所有进程 */
            MPI_Barrier( MPI_COMM_WORLD);
            
            /* 0号进程每隔spt(p-1)个元素选取1个排序后的代表，并将其保存到temp1的[1...spt(p-1)]位置，避开0是因为会选取0位置元素 */
            for (i=1;i<Spt+1;i++)
                temp1[i] = temp1[i*Spt-1];
			/* 处理器P0将这P-1个主元播送到所有处理器中，对应于算法13.5步骤（4）*/
            MPI_Bcast(temp1, sizeof(int)*(1+Spt), MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else
        {
            /* 非0编号进程将temp1中抽样的spt(p-1)代表元素发送给0号进程，使用该进程的ID和偏移ALLTOONE_TYPE作为消息标签(区分别的消息) */
            MPI_Send(temp1,sizeof(int)*Spt,MPI_CHAR,0,ALLTOONE_TYPE+MyID, MPI_COMM_WORLD);
            /*非0编号进程到达该位置，0号进程到达第一个barrier */
            MPI_Barrier( MPI_COMM_WORLD);
            /* 等待0号进程将代表元素排序完成，并同步进程 */
            MPI_Barrier( MPI_COMM_WORLD);
            /* 将选取后的主元广播到所有进程， 并进行同步 */
            MPI_Bcast(temp1, sizeof(int)*(1+Spt), MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }

		/*每个处理器根据上步送来的P-1个主元把自己的n/P个数据分成P段，记为处理器Pi的第j+1段，其中i=0,…,P-1，j=0,…,P-1，对应于算法13.5步骤（5）*/
        n = mylength;
        /* index[2*i]指明段首，index[2*i+1]段末，第一段开始为0 */
        index[0] = 0;
        i = 1;
        /* temp1[i]是分段的界限，此时arr是有序的，即arr[0]是最小的元素，以下循环使arr[0]与各个分段界限进行比较，找到arr[0]属于分段，之前的分段位置都是0 */
        while ((arr[0]>=temp1[i])&&(i<SumID))
        {
            index[2*i-1] = 0;
            index[2*i] = 0;
            i++;
        }
        /* 若arr[0]大于所有分段的元素，则将最后一个位置2p-1设置为n */
        if (i==SumID) index[2*i-1] = n;
        c1 = 0;
        /* 若存在某一个分解元素使得arr[0]小于它，则需要对其进行分段 */
        while (i<SumID)
        {
            /* c4指明当前分界元素，需要在arr找到一个元素等于c4 */
            c4 = temp1[i];
            /* c3为该进程数据长度 */
            c3 = n;
            /* c2为c1和c3的中点，元素是有序的可以二分搜索，查找元素等于c4 */
            c2 = (int)((c1+c3)/2);
            while ((arr[c2]!=c4)&&(c1<c3))
            {
                if (arr[c2]>c4)
                {
                    c3 = c2-1;
                    c2 = (int)((c1+c3)/2);
                }
                else
                {
                    c1 = c2+1;
                    c2 = (int)((c1+c3)/2);
                }
            }
            /* 若找到的元素小等于c4，则将c2指向大于c4的第一个元素 */
            while ((arr[c2]<=c4)&&(c2<n)) c2++;
            /* 若当前所有的元素都小于c4，则将该段末尾分解设置为n，并且将其后的段的首末都设置为0 */
            if (c2==n)
            {
                index[2*i-1] = n;
                for (k=i;k<SumID;k++)
                {
                    index[2*k] = 0;
                    index[2*k+1] = 0;
                }
                i = SumID;
            }
            else
            {
                index[2*i] = c2;
                index[2*i-1] = c2;
            }
            /* c1成为上未分段元素 */
            c1 = c2;
            /* c2取未分段元素中间位置 */
            c2 = (int)((c1+c3)/2);
            /* 进行下一个分界元素检查 */
            i++;
        }
        /* 最后一段的最后一个位置2p-1设置为 n */
        if (i==SumID) index[2*i-1] = n;
        
        /* 等待所有进程分段完成 */
        MPI_Barrier( MPI_COMM_WORLD);

		/*每个处理器送它的第i+1段给处理器Pi，从而使得第i个处理器含有所有处理器的第i段数据(i=0,…,P-1)，，对应于算法13.5步骤（6）*/

        j = 0;
        /* 处理一共 p 个进程 */
        for (i=0;i<SumID;i++)
        {
            if (i==MyID)
            {
                /* Pi进程所做处理 */
                /* 第i段的长度放到temp1[i] */
                temp1[i] = index[2*i+1]-index[2*i];
                for (n=0;n<SumID;n++)
                    /* 将本进程第n段的长度发给n进程 */
                    if (n!=MyID)
                    {
                        k = index[2*n+1]-index[2*n];
                        MPI_Send(&k, sizeof(int), MPI_CHAR, n, MULTI_LEN+MyID,MPI_COMM_WORLD);
                    }
            }
            else
            {
                /* 非Pi进程，则接受第i进程的信息，即每个x进程都收到所有进程的x信息，并存放到temp1[i]位置 */
                MPI_Recv(&temp1[i], sizeof(int), MPI_CHAR, i,MULTI_LEN+i, MPI_COMM_WORLD, &status[j++]);
            }
        }
        
        /* 同步所有进程 */
        MPI_Barrier(MPI_COMM_WORLD);

        j = 0;
        k = 0;
        l = 0;
        
        /* 依次处理每一个进程的数据 */
        for (i=0;i<SumID;i++)
        {
            /* 同步所有进程 */
            MPI_Barrier(MPI_COMM_WORLD);
            
            /* 若处理当前进程的数据，将第i段的数据拷贝到arr1[k] */
            if (i==MyID)
            {
                for (n=index[2*i];n<index[2*i+1];n++)
                    arr1[k++] = arr[n];
            }
            
            /* 同步所有进程 */
            MPI_Barrier(MPI_COMM_WORLD);
            
            /* 若处理当前进程的数据，将第n段的数据发送给n进程 */
            if (i==MyID)
            {
                for (n=0;n<SumID;n++)
                    if (n!=MyID)
                    {
                        MPI_Send(&arr[index[2*n]], sizeof(int)*(index[2*n+1]-index[2*n]),MPI_CHAR, n, MULTI_TYPE+MyID, MPI_COMM_WORLD);
                    }

            }
            else
            {
                l = temp1[i];
                /* 接收来自i进程的数据，并放到arr1[k]位置，接收的是i进程，数据时当前进程需要处理的段 */
                MPI_Recv(&arr1[k], l*sizeof(int), MPI_CHAR, i ,MULTI_TYPE+i, MPI_COMM_WORLD, &status[j++]);
                k=k+l;
            }
            
            /* 同步所有进程 */
            MPI_Barrier(MPI_COMM_WORLD);
        }
        /* 重新设置当前进程的数据长度 */
        mylength = k;
        MPI_Barrier(MPI_COMM_WORLD);

		/*每个处理器再通过P路归并排序将上一步的到的数据排序；从而这n个数据便是有序的，，对应于算法13.5步骤（7） */
        k = 0;
        multimerge(arr1,temp1,arr,&k,SumID);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    printf("On node %d the sorted data is : \n",MyID);
    /* 输出排序后的数据 */
    for (i=0;i<mylength;i++){
        if(k % 2 == 0) printf("%d : ",arr1[i]);
        else printf("%d : ",arr[i]);
    }
    printf("\n");
}


/*输出错误信息*/
merror(char* ch)
{
    printf("%s\n",ch);
    exit(1);
}


/*串行快速排序算法*/
quicksort(int *datas,int bb,int ee)
{
    int tt,i,j;
    tt = datas[bb];
    i = bb;
    j = ee;

    if (i<j)
    {
        while(i<j)
        {
            /* 从右边开始寻找第一个小于tt的位置j */
            while ((i<j)&&(tt<=datas[j])) j--;
            if (i<j)
            {
                /* 将第一个小于tt的元素放到i位置，左侧 */
                datas[i] = datas[j];
                i++;
                /* 从左侧寻找第一个大于等于tt位置i的元素 */
                while ((i<j)&&(tt>datas[i])) i++;
                if (i<j)
                {
                    /* 把i位置的元素放到之前的j */
                    datas[j] = datas[i];
                    j--;
                    /* 若j == i 则将tt放置到i */
                    if (i==j) datas[i] = tt;
                }
                /* 没找到i的位置，则将tt放置到j */
                else datas[j] = tt;
            } else datas[i] = tt;
        }

        quicksort(datas,bb,i-1);
        quicksort(datas,i+1,ee);
    }
}


/*串行多路归并算法*/
multimerge(int *data1,int *ind,int *data,int *iter,int SumID)
{
    int i,j,n;

    j = 0;
    /* 过滤长度为 0 的分段 */
    for (i=0;i<SumID;i++)
        if (ind[i]>0)
        {
            ind[j++] = ind[i];
            if (j<i+1) ind[i] = 0;
        }
    
    /* 若分段数目为两个以上 */
    if ( j>1 )
    {
        n = 0;
        for (i=0;i<j,i+1<j;i=i+2)
        {
            /* 依次合并相邻的两段 */
            merge(&(data1[n]),ind[i],ind[i+1],&(data[n]));
            /* 修改段数对应的长度 */
            ind[i] += ind[i+1];
            ind[i+1] = 0;
            /* 增加n的偏移 */
            n += ind[i];
        }
        /* 多余一段，将其拷贝到data */
        if (j%2==1 )
            for (i=0;i<ind[j-1];i++) data[n++]=data1[n];
        (*iter)++;
        /* 交换data与data1再次进行归并 */
        multimerge(data,ind,data1,iter,SumID);
    }
}


merge(int *data1,int s1,int s2,int *data2)
{
    int i,l,m;

    l = 0;
    m = s1;
    /* l指向一块，m指向另一块 */
    for (i=0;i<s1+s2;i++)
    {
        if (l==s1)
            /* 第一块数据扫描完毕 */
            data2[i]=data1[m++];
        else
            /* 第二块数据扫描完毕 */
            if (m==s2+s1)
                data2[i]=data1[l++];
            else
                /* 选择两块数据中较小的一个元素 */
                if (data1[l]>data1[m])
                    data2[i]=data1[m++];
                else
                    data2[i]=data1[l++];
    }
}
