#include "new_ga.h"
#include "mpi.h"
using namespace std;

int main(int argc,char**argv)
{
	int myid,numprocs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);  //current process id
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); //num of process
	MPI_Status status;
	cout<<numprocs<<endl;
	srand((unsigned int)time(NULL));
	double **sendbuff = NULL;
	double *obj_mpi = NULL;
	double *slaverecv = NULL;
	double *obj_slave = NULL;
	double **sends = NULL;
	double *recvs = NULL;
	double *obj_slaves = NULL;
	double *obj_mpis = NULL;
	//int n=1;
	int times = 900;
	double beta = 0;
	int ll = P/(numprocs -1);
	readData();
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0)
	{
		cout<<"lllllll"<<endl;
		sendbuff = new double*[numprocs];
		sends = new double*[numprocs];
		obj_mpi = new double[ll];
		obj_mpis = new double[ll];
		for(int j=0;j<numprocs;j++)
		{
			sendbuff[j] = new double[ll*PNUM];
			sends[j] = new double[ll*PNUM];
		}
		cout<<"kkkkkkkkkk"<<endl;
		init();
		popinit();
		encoding();
		for(int x = 1;x<numprocs;x++)
		{
			for(int y=0;y<ll;y++)
			{
				for(int k=0;k<PNUM;k++)
				{
					sends[x][y*PNUM+k] = paras[(x-1)*ll+y][k];
				}
			}
		}//for
		for(int x=1;x<numprocs;x++)
		{
			MPI_Send(sends[x],ll*PNUM,MPI_DOUBLE,x,99,MPI_COMM_WORLD);
		}
	}//if
	else
	{
		recvs = new double[ll*PNUM];
		obj_slaves = new double[ll];
		MPI_Recv(recvs,ll*PNUM,MPI_DOUBLE,0,99,MPI_COMM_WORLD,&status);
		for(int i=0;i<ll;i++)
		{
			double temps[PNUM];
			for(int j=0;j<PNUM;j++)
			{
				temps[j] = recvs[i*PNUM+j];
			}
			obj_slaves[i] = objective(temps);
		}
		MPI_Send(obj_slaves,ll,MPI_DOUBLE,0,99,MPI_COMM_WORLD);
	}
	if(myid==0)
	{
		for(int x=1;x<numprocs;x++)
		{
			MPI_Recv(obj_mpis,ll,MPI_DOUBLE,x,99,MPI_COMM_WORLD,&status);
			for(int y=0;y<ll;y++)
			{
				obj[ll*(x-1)+y] = obj_mpis[y];
				//cout<<"obj "<<P+ll*(x-1)+y<<" "<<obj[P+ll*(x-1)+y]<<endl;
			}
		}
	}
	for(int j=0;j<times;j++)
	{
		if(myid==0)
		{
			mutation();
			variation();
			decoding();
			for(int x = 1;x<numprocs;x++)
			{
				for(int y=0;y<ll;y++)
				{
					for(int k=0;k<PNUM;k++)
					{
						sendbuff[x][y*PNUM+k] = paras2[P+(x-1)*ll+y][k];
					}
				}
			}//for
			for(int x=1;x<numprocs;x++)
			{
				MPI_Send(sendbuff[x],ll*PNUM,MPI_DOUBLE,x,99,MPI_COMM_WORLD);
			}
		}//myid==0
		else
		{
			slaverecv = new double[ll*PNUM];
			obj_slave = new double[ll];
			MPI_Recv(slaverecv,ll*PNUM,MPI_DOUBLE,0,99,MPI_COMM_WORLD,&status);
			for(int i=0;i<ll;i++)
			{
				double temps[PNUM];
				for(int j=0;j<PNUM;j++)
				{
					temps[j] = slaverecv[i*PNUM+j];
				}
				obj_slave[i] = objective(temps);
			}
			MPI_Send(obj_slave,ll,MPI_DOUBLE,0,99,MPI_COMM_WORLD);
		}
		if(myid==0)
		{
			for(int x=1;x<numprocs;x++)
			{
				MPI_Recv(obj_mpi,ll,MPI_DOUBLE,x,99,MPI_COMM_WORLD,&status);
				for(int y=0;y<ll;y++)
				{
					obj[P+ll*(x-1)+y] = obj_mpi[y];
				}
			}
			choose();
			cout<<j<<" best = "<<best<<endl;
			fout<<j<<" "<<best<<" ";
			for(int i=0;i<PNUM;i++)
			{
				fout<<paras[0][i]<<" ";
			}
			fout<<endl;
			if(best<=beta)
			{
				break;
			}
		}
	}//times
	fout.close();
	cout<<"best "<<best<<endl;
	MPI_Finalize();
	return 0;
}
