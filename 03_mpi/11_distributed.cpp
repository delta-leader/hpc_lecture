#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

struct Body {
  double x, y, m, fx, fy;
};

// remove force to reduce transfer size
struct Body_transfer {
  double x, y, m;
};

int main(int argc, char** argv) {
  const int N = 20;
  MPI_Init(&argc, &argv);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
<<<<<<< HEAD
  Body ibody[N/size];
  Body_transfer send_body[N/size], recv_body[N/size];
=======
  Body ibody[N/size], jbody[N/size], buffer[N/size];
>>>>>>> f743798ff25f63cf544466b630c34b35525ca76f
  srand48(rank);
  for(int i=0; i<N/size; i++) {
    ibody[i].x = send_body[i].x = drand48();
    ibody[i].y = send_body[i].y = drand48();
    ibody[i].m = send_body[i].m = drand48();
    ibody[i].fx = ibody[i].fy = 0;
  }
  int send_to = (rank - 1 + size) % size;
  MPI_Datatype MPI_BODY;
  MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_BODY);
  MPI_Type_commit(&MPI_BODY);
  // create window
  MPI_Win win;
  // receive data to recv_body
  MPI_Win_create(recv_body, N/size*sizeof(Body_transfer), sizeof(Body_transfer), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  for(int irank=0; irank<size; irank++) {
    // put data into window
    MPI_Win_fence(0, win);
    MPI_Put(send_body, N/size, MPI_BODY, send_to, 0, N/size, MPI_BODY, win);
    MPI_Win_fence(0, win);
    for(int i=0; i<N/size; i++) {
      // calculate forces with received data
      for(int j=0; j<N/size; j++) {
        double rx = ibody[i].x - recv_body[j].x;
        double ry = ibody[i].y - recv_body[j].y;
        double r = std::sqrt(rx * rx + ry * ry);
        if (r > 1e-15) {
          ibody[i].fx -= rx * recv_body[j].m / (r * r * r);
          ibody[i].fy -= ry * recv_body[j].m / (r * r * r);
        }
      }
      // copy processed data to send_body
      send_body[i].x = recv_body[i].x;
      send_body[i].y = recv_body[i].y;
      send_body[i].m = recv_body[i].m;
    }
  }
  // free window
  MPI_Win_free(&win);
  for(int irank=0; irank<size; irank++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if(irank==rank) {
      for(int i=0; i<N/size; i++) {
        printf("%d %g %g\n",i+rank*N/size,ibody[i].fx,ibody[i].fy);
      }
    }
  }
  MPI_Finalize();
}
