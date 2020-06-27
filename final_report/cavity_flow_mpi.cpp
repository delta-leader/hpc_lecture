#include <iostream>
#include "mpi.h"

void print(double* data, int ny, int nx){
for(int i=0; i<ny; ++i){
    std::cout<<i<<": [";
    for(int j=0; j<nx; ++j)
        std::cout<<data[i*nx+j]<<" ";
    std::cout<<"]"<<std::endl;
    }
}

void print(double * data, int ny, int nx, int rank, int size){
    for(int irank=0; irank<size; irank++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if(irank==rank) {
	    int offset = rank==0?0:1;
	    for(int i=offset; i<((rank==size-1)?ny:ny-1); ++i){
    		std::cout<<i+rank*size-offset<<"("<<rank<<")"<<": [";
    		for(int j=0; j<nx; ++j)
        		std::cout<<data[i*nx+j]<<" ";
    			std::cout<<"]"<<std::endl;
    		}
      }
    }
}

void print_all(double *data, int ny, int nx, int rank, int size){
    for(int irank=0; irank<size; irank++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if(irank==rank) {
	    for(int i=0; i<ny; ++i){
    		std::cout<<i+rank*size-(rank==0?0:1)<<"("<<rank<<")"<<": [";
    		for(int j=0; j<nx; ++j)
        	    std::cout<<data[i*nx+j]<<" ";
    		std::cout<<"]"<<std::endl;
    	    }
        }
    }
}

void communicate_data(double *sendbuf, double *recvbuf, int ny, int nx, int rank, int size){
    //we never use the data stored in the first and last column, so no need to send it
    int data_count = nx - 2;
    //offsets for sending and receiving - ommits first column
    int sendFirstRow = nx + 1;
    int recvFirstRow = 1;
    int sendLastRow = (ny-2)*nx + 1;
    int recvLastRow = (ny-1)*nx + 1;

    //target processes
    int firstRow = rank - 1;
    int lastRow = rank + 1;

    // the first process only exchanges his last row
    if (rank == 0){
        MPI_Send(sendbuf+sendLastRow, data_count, MPI_DOUBLE, lastRow, 0, MPI_COMM_WORLD);
        MPI_Recv(recvbuf+recvLastRow, data_count, MPI_DOUBLE, lastRow, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
   
    // processes in the middle exchange first & last rows	
    if ((rank > 0) && (rank < size-1)){
        MPI_Sendrecv(sendbuf+sendLastRow, data_count, MPI_DOUBLE, lastRow, 0, recvbuf+recvFirstRow, data_count, MPI_DOUBLE, firstRow, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendbuf+sendFirstRow, data_count, MPI_DOUBLE, firstRow, 0, recvbuf+recvLastRow, data_count, MPI_DOUBLE, lastRow, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // the last process only exhanges the first row
    if (rank == size-1){
        MPI_Recv(recvbuf+recvFirstRow, data_count, MPI_DOUBLE, firstRow, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(sendbuf+sendFirstRow, data_count, MPI_DOUBLE, firstRow, 0, MPI_COMM_WORLD);
    }
}

void communicate_data(double *buf, int ny, int nx, int rank, int size){
    communicate_data(buf, buf, ny, nx, rank, size);
}

void build_up_b(double* b, double rho, double dt, double* u, double* v, double dx, double dy, int nx, int ny){
    for(int i=1; i<ny-1; ++i){
        for(int j=1; j<nx-1; ++j){
            b[i*nx+j] = rho * (1/dt * ((u[i*nx+j+1]-u[i*nx+j-1])  / (2*dx) + (v[(i+1)*nx+j]-v[(i-1)*nx+j]) / (2*dy))
                               -((u[i*nx+j+1]-u[i*nx+j-1]) / (2*dx)) * ((u[i*nx+j+1]-u[i*nx+j-1]) / (2*dx))
                               -2 * (u[(i+1)*nx+j]-u[(i-1)*nx+j]) / (2*dy) * (v[i*nx+j+1]-v[i*nx+j-1]) / (2*dx)
                               -((v[(i+1)*nx+j]-v[(i-1)*nx+j]) / (2*dy)) * ((v[(i+1)*nx+j]-v[(i-1)*nx+j]) / (2*dy))
                              );
        }
           
    }
}

void pressure_poisson_parallel(double* p, double* b, double dx, double dy, int nx, int ny, int nit, int rank, int size){

    double pn[ny*nx];
    for (int t=0; t<nit; ++t){
        //copy old values
        for(int i=0; i<nx*ny; ++i)
            pn[i] = p[i];
        for(int i=1; i<ny-1; ++i)
            for(int j=1; j<nx-1; ++j){
                p[i*nx+j] = ((pn[i*nx+j+1]+pn[i*nx+j-1])*dy*dy + (pn[(i+1)*nx+j]+pn[(i-1)*nx+j])*dx*dx)
                             / (2 * (dx*dx + dy*dy))
                             - dx*dx*dy*dy / (2 * (dx*dx + dy*dy)) * b[i*nx+j];
            }

        //boundary conditions for columns
        for(int i=0; i<ny; ++i){
            p[i*nx+nx-1] = p[i*nx+nx-2];
            p[i*nx] = p[i*nx+1];
        }
        //boundary conditions for rows only for first and last process
        if (rank == 0)
            for (int j=0; j<nx; ++j)
                p[j] = p[nx+j];

	//this can be ommitted because we initalize the last row to 0 and never change it
	//it's just here for completeness sake
        /*if (rank == size-1)
            for (int j =0; j<nx; ++j)
                p[(ny-1)*nx+j] = 0;
	*/

	// update p across processes
        communicate_data(p, ny, nx, rank, size);
    }   
} 

void cavity_flow_parallel(int nt, double* u, double* v, double dt, double dx, double dy, double* p, double rho, double nu, double* b, int nx, int ny, int nit, int rank, int size){

    //setting the row boundary conditions (initially 0)
    //moving this outside of the loop gets rid of the first wasteful round of computations were everything is 0 anyway.
    if (rank == size-1)
        for(int j=0; j<nx; ++j)
            u[(ny-1)*nx+j] = 1;

    // before we start we make sure that all values of p are updated
    // this is not strictly necessary because it is initalized to 0 in the beginning, but nice to have if the initalization changes
    communicate_data(p, ny, nx, rank, size);
    // we don't need to do this for u and v because they get updated at the beginning of each iteration (un, vn)    

    double un[ny*nx], vn[ny*nx];
    for(int n=0; n<nt; ++n){
        //copy old values
        for(int i=0; i<ny*nx; ++i){
            un[i] = u[i];
            vn[i] = v[i];
        }
        // directly update the buffers used for updating u and v (because we don't care about the values of the sent rows in the final result
	communicate_data(u, un, ny, nx, rank, size);
	communicate_data(v, vn, ny, nx, rank, size);
        
	build_up_b(b, rho, dt, un, vn, dx, dy, nx, ny);
        pressure_poisson_parallel(p, b, dx, dy, nx, ny, nit, rank, size);

	//update u and v
	for(int i=1; i<ny-1; ++i)
            for(int j=1; j<nx-1; ++j){
                u[i*nx+j] = un[i*nx+j]-un[i*nx+j]* dt/dx * (un[i*nx+j] - un[i*nx+j-1])
                            - vn[i*nx+j] * dt/dy * (un[i*nx+j] - un[(i-1)*nx+j])
                            - dt / (2*rho*dx) * (p[i*nx+j+1] - p[i*nx+j-1])
                            + nu * (dt/(dx*dx) * (un[i*nx+j+1] - 2*un[i*nx+j] + un[i*nx+j-1])
                            + dt / (dy*dy) * (un[(i+1)*nx+j] - 2*un[i*nx+j] + un[(i-1)*nx+j]));

                v[i*nx+j] = vn[i*nx+j]-un[i*nx+j]* dt/dx * (vn[i*nx+j] - vn[i*nx+j-1])
                            - vn[i*nx+j] * dt/dy * (vn[i*nx+j] - vn[(i-1)*nx+j])
                            - dt / (2*rho*dy) * (p[(i+1)*nx+j] - p[(i-1)*nx+j])
                            + nu * (dt/(dx*dx) * (vn[i*nx+j+1] - 2*vn[i*nx+j] + vn[i*nx+j-1])
                            + dt / (dy*dy) * (vn[(i+1)*nx+j] - 2*vn[i*nx+j] + vn[(i-1)*nx+j]));
        }

	//no need to enforce column boundary conditions (0-initalized)
        //here for completeness sake
        /*for(int i=0; i<ny; ++i){
            u[i*nx] = u[i*nx+nx-1] = 0;
            v[i*nx] = v[i*nx+nx-1] = 0;
        }*/
	//row boundary conditions do not need to be enforced every iteration because those rows are never updated. Just setting them in the beginning is enough	
	//here for completness sake
        /*if (rank == 0)
	    for(int j=0; j<nx; ++j){
                u[j] = 0;
                v[j] = 0;
        }
	if (rank == size-1)
            for(int j=0; j<nx; ++j){
                u[(ny-1)*nx+j] = 1;
                v[(ny-1)*nx+j] = 0;
        }*/
    }
}

void pressure_poisson_serial(double* p, double* b, double dx, double dy, int nx, int ny, int nit){
    double pn[ny*nx];
    for (int t=0; t<nit; ++t){
        //copy old values
        for(int i=0; i<nx*ny; ++i)
            pn[i] = p[i];

        for(int i=0; i<ny-1; ++i)
            for(int j=0; j<nx-1; ++j){
                p[i*nx+j] = ((pn[i*nx+j+1]+pn[i*nx+j-1])*dy*dy + (pn[(i+1)*nx+j]+pn[(i-1)*nx+j])*dx*dx)
                             / (2 * (dx*dx + dy*dy))
                             - dx*dx*dy*dy / (2 * (dx*dx + dy*dy)) * b[i*nx+j];
            }
        
        //boundary conditions
        for(int i=0; i<ny; ++i){
            p[i*nx+nx-1] = p[i*nx+nx-2];
            p[i*nx] = p[i*nx+1];
        }
        for (int j =0; j<nx; ++j){
            p[j] = p[nx+j];
            // not necessary because 0 initialized
            //p[(ny-1)*nx+j] = 0; 
        }
    } 
} 

void cavity_flow_serial(int nt, double* u, double* v, double dt, double dx, double dy, double* p, double rho, double nu, double* b, int nx, int ny, int nit){

    //setting the boundary condition in the beginning gets rid of one round of wasteful computations
    for(int j=0; j<nx; ++j)
        u[(ny-1)*nx+j] = 1;

    double un[ny*nx], vn[ny*nx];
    for(int n=0; n<nt; ++n){

        //copy old values
        for(int i=0; i<ny*nx; ++i){
            un[i] = u[i];
            vn[i] = v[i];
        }
        build_up_b(b, rho, dt, u, v, dx, dy, nx, ny);
        pressure_poisson_serial(p, b, dx, dy, nx, ny, nit);

        for(int i=1; i<ny-1; ++i)
            for(int j=1; j<nx-1; ++j){
                u[i*nx+j] = un[i*nx+j]-un[i*nx+j]* dt/dx * (un[i*nx+j] - un[i*nx+j-1])
                            - vn[i*nx+j] * dt/dy * (un[i*nx+j] - un[(i-1)*nx+j])
                            - dt / (2*rho*dx) * (p[i*nx+j+1] - p[i*nx+j-1])
                            + nu * (dt/(dx*dx) * (un[i*nx+j+1] - 2*un[i*nx+j] + un[i*nx+j-1])
                            + dt / (dy*dy) * (un[(i+1)*nx+j] - 2*un[i*nx+j] + un[(i-1)*nx+j]));

                v[i*nx+j] = vn[i*nx+j]-un[i*nx+j]* dt/dx * (vn[i*nx+j] - vn[i*nx+j-1])
                            - vn[i*nx+j] * dt/dy * (vn[i*nx+j] - vn[(i-1)*nx+j])
                            - dt / (2*rho*dy) * (p[(i+1)*nx+j] - p[(i-1)*nx+j])
                            + nu * (dt/(dx*dx) * (vn[i*nx+j+1] - 2*vn[i*nx+j] + vn[i*nx+j-1])
                            + dt / (dy*dy) * (vn[(i+1)*nx+j] - 2*vn[i*nx+j] + vn[(i-1)*nx+j]));
        }
        //because of initalization we can ommit boundary conditions
        //boundary conditions
        /*for(int i=0; i<ny; ++i){
            u[i*nx] = u[i*nx+nx-1] = 0;
            v[i*nx] = v[i*nx+nx-1] = 0;
        }   
        for(int j=0; j<nx; ++j){
            u[j] = 0;
            u[(ny-1)*nx+j] = 1;
            v[j] = v[(ny-1)*nx+j] = 0;
        }*/
    }
}

int main(int argc, char** argv) {

    //settings
    int nx = 41;
    int ny = 41;
    //one iteration less than in the python example because I moved the boundary condition out of the loop
    int nt = 99;
    int nit = 50;
    int c = 1;
    double dx = 2 / (nx - 1.);
    double dy = 2 / (ny - 1.);
    double rho = 1;
    double nu = 0.1;
    double dt = 0.001;

    //initialize mpi
    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = ny/size;
    int remainder = ny%size;

    //only start a distributed execution if there are at least 2 rows perprocess
    //otherwise the communication overhead is too much (2 communications per processed row)
    if (count > 2){
        int rows = (rank < remainder) ? count+1 : count;
        int comm_rows = (rank==0 || rank==size-1) ? 1 : 2;
        int array_size = (rows + comm_rows) * nx;
        //printf("rank: %d/%d, %d-%d,%d/%d\n",rank,size,start, stop, rows, array_size);	

        double u[array_size], v[array_size], p[array_size], b[array_size];
        
        //initalize to 0
        for (int i=0; i<array_size; ++i)
            u[i] = v[i] = p[i] = b[i] = 0;

        //start distributed calculation
        cavity_flow_parallel(nt, u, v, dt, dx, dy, p, rho, nu, b, nx, rows+comm_rows, nit, rank, size);

        //print distributed data
        //print_all(p, rows+comm_row, nx, rank, size);
        //print(u, rows+comm_row, nx, rank, size);

        //gather the data on the root process
        double u_all[nx*ny], v_all[nx*ny], p_all[nx*ny];
        if (rank == 0){
	    int displs[size], rcounts[size];
	    int offset =0;
	    for (int i=0; i<size; ++i){
                rcounts[i] = (i<remainder?count+1:count) * nx;
	        displs[i] = offset;
	        offset += rcounts[i];
            }
            MPI_Gatherv(u, rows*nx, MPI_DOUBLE, u_all, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv(v, rows*nx, MPI_DOUBLE, v_all, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv(p, rows*nx, MPI_DOUBLE, p_all, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        else {
            MPI_Gatherv(u+nx, rows*nx, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv(v+nx, rows*nx, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv(p+nx, rows*nx, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

	//print gathered data
        if (rank==0)
	    print(v_all, ny, nx);
    }

    else {
        if (rank==0){
            printf("Parallel overhead was too large, computation was executed in serial.\n");

            double u[ny*nx], v[ny*nx], p[ny*nx], b[ny*nx];        
            //initalize to 0
            for (int i=0; i<ny*nx; ++i)
                u[i] = v[i] = p[i] = b[i] = 0;

            cavity_flow_serial(nt, u, v, dt, dx, dy, p, rho, nu,  b, nx, ny, nit);
	    print(u, ny, nx);
        }
    }

    MPI_Finalize();

}
