#ifndef _H_CLUSTER_H_
#define _H_CLUSTER_H_

extern int cluster[50000];
extern int cluster_atoms[50000];
extern int total_cluster_atoms;
extern int noncluster_atoms[50000];
extern int total_noncluster_atoms;

//size of the cluster
extern int cluster_size[1000];

//neighbors of the current atom
extern int curr_neigh[5];

// coincedence flaq of the current column
extern int coin_flaq;

extern int total_atom_number;

//cluster matrix
extern int cluster_direct[5000][5000];
extern int cluster_sum[5000][5000];

extern int cluster_column[20010];
extern int cluster_function[10000];
extern int curr_size [10000];
//coordinates of the atoms in the system
extern double sys_adcoords[90000][3];
extern double atoms_adcoords[90000][3];

extern int row_index;
extern int row_coin_flaq;

extern FILE *F_cluster;

void DataRead( void );

void MatrixMade( void );

int RowCompare( int a, int b );

void ClusterData( void );

void ClusterDistribution( void );

//void ColumnEqual( void );


#endif
