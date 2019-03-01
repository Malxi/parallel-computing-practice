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
    /* ��ȡ��ǰ���̱�� */
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    
    /* PlusNum Ϊ60,���ݴ�СΪBaseNum*PlusNum */
    PlusNum=60;
    DataSize = BaseNum*PlusNum;
    
    /* ��Ž������ PlusNum */
    if (MyID==0)
        printf("The DataSize is : %lu\n",PlusNum);
    
    /* Psrs_Main���� */
    Psrs_Main();
    
    /* ��Ž������\n */
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
    
    /* MyId�ǵ�ǰ���̱��,SumID�ǽ�����Ŀ */
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);
    MPI_Comm_size(MPI_COMM_WORLD,&SumID);
    
    /* SptΪ���һ�����̱�� */
    Spt = SumID-1;

	/*��ʼ������*/
    /* ��������DataSize��С�ڴ�,arrָ��ʼ,arr1ָ��ڶ��ݿ�ʼ */
    arr = (int *)malloc(2*DataSize*sizeof(int));
    if (arr==0) merror("malloc memory for arr error!");
    arr1 = &arr[DataSize];
    
    /* ������������1,�����SumID*Spt�ڴ��temp1,2*SumID��index, temp1�洢����Ԫ�ش�С����ΪSumID*Spt�� index�����δ����СΪ2*SumID */
    if (SumID>1)
    {
        temp1 = (int *)malloc(sizeof(int)*SumID*Spt);
        if (temp1==0) merror("malloc memory for temp1 error!");
        index = (int *)malloc(sizeof(int)*2*SumID);
        if (index==0) merror("malloc memory for index error!");
    }
    
    /* �ȴ����н��̶���������ʼ��,�ڴ������� */
    MPI_Barrier( MPI_COMM_WORLD);
    
    /* ���㵱ǰ�������򳤶�,�������������Ϊ��ǰ���̱�� */
    mylength = DataSize / SumID;
    srand(MyID);
    
    /* �����ǰ�ڵ㼰�ýڵ����������arr,�����������,���������ǹ̶���,�ڽڵ�����ͬʱ,��������ͬ�� */
    printf("This is node %d \n",MyID);
    printf("On node %d the input data is:\n",MyID);
    for (i=0;i<mylength;i++)
    {
        arr[i] = (int)rand();
        printf("%d : ",arr[i]);
    }
    printf("\n");

	/*ÿ�����������Լ���n/P�������ô��п�������(Quicksort)���õ�һ���ź�������У���Ӧ���㷨13.5���裨1��*/
    MPI_Barrier( MPI_COMM_WORLD);
    quicksort(arr,0,mylength - 1);
    MPI_Barrier( MPI_COMM_WORLD);

	/*ÿ�����������ź����������ѡȡ��w��2w��3w������(P-1)w����P-1��������Ϊ����Ԫ�أ�����w=n/P*P����Ӧ���㷨13.5���裨2��*/
    if (SumID>1)
    {
        /* �ȴ����н��̶����Լ�������������� */
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* nΪ�����ļ����ÿ�����̶���ȡspt(p-1)������Ԫ�أ�������ŵ�temp1 */
        n = (int)(mylength/(Spt+1));
        for (i=0;i<Spt;i++)
            temp1[i] = arr[(i+1)*n-1];

        /* �ȴ����н��̳������ */
        MPI_Barrier(MPI_COMM_WORLD);

        if (MyID==0)
        {
			/* ÿ����������ѡ�õĴ���Ԫ���͵�������P0�У���Ӧ���㷨13.5���裨3�� */
            /* Irecv����Ҫ�ȴ����ͻ������Ϣ��ɾͿ���ִ���������񣬲������õ�request�������У�ÿ������һ���ֲ� */
            /* �������ݷ�����temp1��ÿ��sptƫ��λ�� */
            j = 0;
            for (i=1;i<SumID;i++)
                MPI_Irecv(&temp1[i*Spt], sizeof(int)*Spt, MPI_CHAR, i,ALLTOONE_TYPE+i, MPI_COMM_WORLD, &request[j++]);
            /* �ȴ�����������̵Ĵ���Ԫ�� */
            MPI_Waitall(SumID-1, request, status);

			/* ������P0����һ��������P�����������������P·�鲢����ѡ�������ĵ�P-1��2(P-1)������(P-1)(P-1)����P-1����Ԫ������Ӧ���㷨13.5���裨3��*/
            /* �ȴ�������̵����λ�� */
            MPI_Barrier(MPI_COMM_WORLD);
            quicksort(temp1,0,SumID*Spt-1);
            /* ͬ�����н��� */
            MPI_Barrier( MPI_COMM_WORLD);
            
            /* 0�Ž���ÿ��spt(p-1)��Ԫ��ѡȡ1�������Ĵ��������䱣�浽temp1��[1...spt(p-1)]λ�ã��ܿ�0����Ϊ��ѡȡ0λ��Ԫ�� */
            for (i=1;i<Spt+1;i++)
                temp1[i] = temp1[i*Spt-1];
			/* ������P0����P-1����Ԫ���͵����д������У���Ӧ���㷨13.5���裨4��*/
            MPI_Bcast(temp1, sizeof(int)*(1+Spt), MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else
        {
            /* ��0��Ž��̽�temp1�г�����spt(p-1)����Ԫ�ط��͸�0�Ž��̣�ʹ�øý��̵�ID��ƫ��ALLTOONE_TYPE��Ϊ��Ϣ��ǩ(���ֱ����Ϣ) */
            MPI_Send(temp1,sizeof(int)*Spt,MPI_CHAR,0,ALLTOONE_TYPE+MyID, MPI_COMM_WORLD);
            /*��0��Ž��̵����λ�ã�0�Ž��̵����һ��barrier */
            MPI_Barrier( MPI_COMM_WORLD);
            /* �ȴ�0�Ž��̽�����Ԫ��������ɣ���ͬ������ */
            MPI_Barrier( MPI_COMM_WORLD);
            /* ��ѡȡ�����Ԫ�㲥�����н��̣� ������ͬ�� */
            MPI_Bcast(temp1, sizeof(int)*(1+Spt), MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }

		/*ÿ�������������ϲ�������P-1����Ԫ���Լ���n/P�����ݷֳ�P�Σ���Ϊ������Pi�ĵ�j+1�Σ�����i=0,��,P-1��j=0,��,P-1����Ӧ���㷨13.5���裨5��*/
        n = mylength;
        /* index[2*i]ָ�����ף�index[2*i+1]��ĩ����һ�ο�ʼΪ0 */
        index[0] = 0;
        i = 1;
        /* temp1[i]�ǷֶεĽ��ޣ���ʱarr������ģ���arr[0]����С��Ԫ�أ�����ѭ��ʹarr[0]������ֶν��޽��бȽϣ��ҵ�arr[0]���ڷֶΣ�֮ǰ�ķֶ�λ�ö���0 */
        while ((arr[0]>=temp1[i])&&(i<SumID))
        {
            index[2*i-1] = 0;
            index[2*i] = 0;
            i++;
        }
        /* ��arr[0]�������зֶε�Ԫ�أ������һ��λ��2p-1����Ϊn */
        if (i==SumID) index[2*i-1] = n;
        c1 = 0;
        /* ������ĳһ���ֽ�Ԫ��ʹ��arr[0]С����������Ҫ������зֶ� */
        while (i<SumID)
        {
            /* c4ָ����ǰ�ֽ�Ԫ�أ���Ҫ��arr�ҵ�һ��Ԫ�ص���c4 */
            c4 = temp1[i];
            /* c3Ϊ�ý������ݳ��� */
            c3 = n;
            /* c2Ϊc1��c3���е㣬Ԫ��������Ŀ��Զ�������������Ԫ�ص���c4 */
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
            /* ���ҵ���Ԫ��С����c4����c2ָ�����c4�ĵ�һ��Ԫ�� */
            while ((arr[c2]<=c4)&&(c2<n)) c2++;
            /* ����ǰ���е�Ԫ�ض�С��c4���򽫸ö�ĩβ�ֽ�����Ϊn�����ҽ����Ķε���ĩ������Ϊ0 */
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
            /* c1��Ϊ��δ�ֶ�Ԫ�� */
            c1 = c2;
            /* c2ȡδ�ֶ�Ԫ���м�λ�� */
            c2 = (int)((c1+c3)/2);
            /* ������һ���ֽ�Ԫ�ؼ�� */
            i++;
        }
        /* ���һ�ε����һ��λ��2p-1����Ϊ n */
        if (i==SumID) index[2*i-1] = n;
        
        /* �ȴ����н��̷ֶ���� */
        MPI_Barrier( MPI_COMM_WORLD);

		/*ÿ�������������ĵ�i+1�θ�������Pi���Ӷ�ʹ�õ�i���������������д������ĵ�i������(i=0,��,P-1)������Ӧ���㷨13.5���裨6��*/

        j = 0;
        /* ����һ�� p ������ */
        for (i=0;i<SumID;i++)
        {
            if (i==MyID)
            {
                /* Pi������������ */
                /* ��i�εĳ��ȷŵ�temp1[i] */
                temp1[i] = index[2*i+1]-index[2*i];
                for (n=0;n<SumID;n++)
                    /* �������̵�n�εĳ��ȷ���n���� */
                    if (n!=MyID)
                    {
                        k = index[2*n+1]-index[2*n];
                        MPI_Send(&k, sizeof(int), MPI_CHAR, n, MULTI_LEN+MyID,MPI_COMM_WORLD);
                    }
            }
            else
            {
                /* ��Pi���̣�����ܵ�i���̵���Ϣ����ÿ��x���̶��յ����н��̵�x��Ϣ������ŵ�temp1[i]λ�� */
                MPI_Recv(&temp1[i], sizeof(int), MPI_CHAR, i,MULTI_LEN+i, MPI_COMM_WORLD, &status[j++]);
            }
        }
        
        /* ͬ�����н��� */
        MPI_Barrier(MPI_COMM_WORLD);

        j = 0;
        k = 0;
        l = 0;
        
        /* ���δ���ÿһ�����̵����� */
        for (i=0;i<SumID;i++)
        {
            /* ͬ�����н��� */
            MPI_Barrier(MPI_COMM_WORLD);
            
            /* ������ǰ���̵����ݣ�����i�ε����ݿ�����arr1[k] */
            if (i==MyID)
            {
                for (n=index[2*i];n<index[2*i+1];n++)
                    arr1[k++] = arr[n];
            }
            
            /* ͬ�����н��� */
            MPI_Barrier(MPI_COMM_WORLD);
            
            /* ������ǰ���̵����ݣ�����n�ε����ݷ��͸�n���� */
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
                /* ��������i���̵����ݣ����ŵ�arr1[k]λ�ã����յ���i���̣�����ʱ��ǰ������Ҫ����Ķ� */
                MPI_Recv(&arr1[k], l*sizeof(int), MPI_CHAR, i ,MULTI_TYPE+i, MPI_COMM_WORLD, &status[j++]);
                k=k+l;
            }
            
            /* ͬ�����н��� */
            MPI_Barrier(MPI_COMM_WORLD);
        }
        /* �������õ�ǰ���̵����ݳ��� */
        mylength = k;
        MPI_Barrier(MPI_COMM_WORLD);

		/*ÿ����������ͨ��P·�鲢������һ���ĵ����������򣻴Ӷ���n�����ݱ�������ģ�����Ӧ���㷨13.5���裨7�� */
        k = 0;
        multimerge(arr1,temp1,arr,&k,SumID);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    printf("On node %d the sorted data is : \n",MyID);
    /* ������������� */
    for (i=0;i<mylength;i++){
        if(k % 2 == 0) printf("%d : ",arr1[i]);
        else printf("%d : ",arr[i]);
    }
    printf("\n");
}


/*���������Ϣ*/
merror(char* ch)
{
    printf("%s\n",ch);
    exit(1);
}


/*���п��������㷨*/
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
            /* ���ұ߿�ʼѰ�ҵ�һ��С��tt��λ��j */
            while ((i<j)&&(tt<=datas[j])) j--;
            if (i<j)
            {
                /* ����һ��С��tt��Ԫ�طŵ�iλ�ã���� */
                datas[i] = datas[j];
                i++;
                /* �����Ѱ�ҵ�һ�����ڵ���ttλ��i��Ԫ�� */
                while ((i<j)&&(tt>datas[i])) i++;
                if (i<j)
                {
                    /* ��iλ�õ�Ԫ�طŵ�֮ǰ��j */
                    datas[j] = datas[i];
                    j--;
                    /* ��j == i ��tt���õ�i */
                    if (i==j) datas[i] = tt;
                }
                /* û�ҵ�i��λ�ã���tt���õ�j */
                else datas[j] = tt;
            } else datas[i] = tt;
        }

        quicksort(datas,bb,i-1);
        quicksort(datas,i+1,ee);
    }
}


/*���ж�·�鲢�㷨*/
multimerge(int *data1,int *ind,int *data,int *iter,int SumID)
{
    int i,j,n;

    j = 0;
    /* ���˳���Ϊ 0 �ķֶ� */
    for (i=0;i<SumID;i++)
        if (ind[i]>0)
        {
            ind[j++] = ind[i];
            if (j<i+1) ind[i] = 0;
        }
    
    /* ���ֶ���ĿΪ�������� */
    if ( j>1 )
    {
        n = 0;
        for (i=0;i<j,i+1<j;i=i+2)
        {
            /* ���κϲ����ڵ����� */
            merge(&(data1[n]),ind[i],ind[i+1],&(data[n]));
            /* �޸Ķ�����Ӧ�ĳ��� */
            ind[i] += ind[i+1];
            ind[i+1] = 0;
            /* ����n��ƫ�� */
            n += ind[i];
        }
        /* ����һ�Σ����俽����data */
        if (j%2==1 )
            for (i=0;i<ind[j-1];i++) data[n++]=data1[n];
        (*iter)++;
        /* ����data��data1�ٴν��й鲢 */
        multimerge(data,ind,data1,iter,SumID);
    }
}


merge(int *data1,int s1,int s2,int *data2)
{
    int i,l,m;

    l = 0;
    m = s1;
    /* lָ��һ�飬mָ����һ�� */
    for (i=0;i<s1+s2;i++)
    {
        if (l==s1)
            /* ��һ������ɨ����� */
            data2[i]=data1[m++];
        else
            /* �ڶ�������ɨ����� */
            if (m==s2+s1)
                data2[i]=data1[l++];
            else
                /* ѡ�����������н�С��һ��Ԫ�� */
                if (data1[l]>data1[m])
                    data2[i]=data1[m++];
                else
                    data2[i]=data1[l++];
    }
}
