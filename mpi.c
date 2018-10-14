#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <fcntl.h>
 
enum mode { // ENUM for different image modes
    GRAY,
    COLOR
};
 
// Convoluton function
// Applies filter on given pixel array input, into the output array
// The convolution is applied on every image pixel in the rectangle created
// by 4 points, startX, startY, endX, endY
void convolution( unsigned char * input, unsigned char * output, double * filter, int rowLength, int startX, int endX, int startY, int endY, enum mode color ) {
    unsigned char * currentPixel = NULL;
    for ( int i = startX; i <= endX; i++ ) {
        for ( int j = startY; j <= endY; j++ ) {
            // For each pixel inside our rectangle, apply convolution
            int filterIndex = 0;
            if ( color == GRAY ) {
                double result = 0;
                for ( int offY = -1; offY <= 1; offY++ ) {
                    for ( int offX = -1; offX <= 1; offX++ ) {
                        // Sum neighbour pixels into result and multiply each
                        // of them with the appropriate filter value
                        currentPixel = &input[ ( i + offX ) + rowLength * ( j + offY ) ];
                        result += *currentPixel * filter[ filterIndex++ ];
                    }
                }
                // Write pixel into the output array
                output[ j * rowLength + i] = result;
            } else if ( color == COLOR ) {
                double red = 0, green = 0 , blue = 0;
                for ( int offY = -1; offY <= 1; offY++ ) {
                    for ( int offX = -1; offX <= 1; offX++ ) {
                        // Sum neighbour pixels into result and multiply each
                        // of them with the appropriate filter value
                        currentPixel = &input[ 3 * ( i + offX ) + 3*rowLength * ( j + offY ) ];
                        red += *currentPixel * filter[ filterIndex ];
                        green += *(currentPixel + sizeof(unsigned char) ) * filter[ filterIndex ];
                        blue += *(currentPixel + 2*sizeof(unsigned char)) * filter[ filterIndex++ ];
                    }
                }
                // Write pixel into the output array
                output[ j * 3*rowLength + 3*i] = red;
                output[ j * 3*rowLength + 3*i + sizeof(unsigned char) ] = green;
                output[ j * 3*rowLength + 3*i + 2*sizeof(unsigned char) ] = blue;
            }
           
 
        }
    }
}
 
 
int main(int argc, char* argv[] ) {
 
    int processes;
    int myRank;
 
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
 
    // Select color mode, store it in <enum color>
    enum mode color = GRAY;
   
    // initialize filters here
    double filter[9] = {0,0,0,0,1,0,0,0,0};    
    //get arguments
    //execute form : mpi width height filename repeats mode
    if (argc != 6){
        perror("Arguments");
        return(0);
        //error
        //exit
    }
    int width = atoi(argv[1]), height = atoi(argv[2]);
    char * fileName;
    if ((fileName = malloc( (strlen(argv[3])+1) *sizeof(char))) == NULL) {
        perror("Malloc");
        return 0;
    }
    strcpy(fileName, argv[3]);
    int repeats  = atoi(argv[4]);
   
    if (!(strcmp(argv[5], "GRAY"))){
        color = GRAY;    
    }
    if (!(strcmp(argv[5], "COLOR"))){
        color = COLOR;    
    }
 
    //initialize neighbours
    int N, S, E = MPI_PROC_NULL, W = MPI_PROC_NULL, NE = MPI_PROC_NULL, NW = MPI_PROC_NULL, SE = MPI_PROC_NULL, SW = MPI_PROC_NULL;
    int root = (int) sqrt(processes);
   
    //send requests
    MPI_Request sendN;
    MPI_Request sendNE;
    MPI_Request sendE;
    MPI_Request sendSE;
    MPI_Request sendS;
    MPI_Request sendSW;
    MPI_Request sendW;
    MPI_Request sendNW;
 
    //recieve reequests
    MPI_Request recN;
    MPI_Request recNE;
    MPI_Request recE;
    MPI_Request recSE;
    MPI_Request recS;
    MPI_Request recSW;
    MPI_Request recW;
    MPI_Request recNW;
   
    //row data type
    MPI_Datatype grayRow, colorRow;
    MPI_Type_contiguous(width/root, MPI_BYTE, &grayRow);
    MPI_Type_commit(&grayRow);
    MPI_Type_contiguous(3*(width/root), MPI_BYTE, &colorRow);
    MPI_Type_commit(&colorRow);
 
    //colum data type
    MPI_Datatype grayCol, colorCol;
    MPI_Type_vector(height/root, 1, (width/root)+2, MPI_BYTE, &grayCol);
    MPI_Type_commit(&grayCol);
    MPI_Type_vector(height/root, 3, 3*((width/root)+2), MPI_BYTE, &colorCol);
    MPI_Type_commit(&colorCol);
 
    //corner data type
    MPI_Datatype grayCorner, colorCorner;
    MPI_Type_contiguous(1, MPI_BYTE, &grayCorner);
    MPI_Type_commit(&grayCorner);
    MPI_Type_contiguous(3, MPI_BYTE, &colorCorner);
    MPI_Type_commit(&colorCorner);
   
    //find neighbours
    //N
    if( (N = myRank - root) < 0){
        N = MPI_PROC_NULL;
    }
   
    //S
    if( (S = myRank + root) >= processes){
        S = MPI_PROC_NULL;
    }
 
    //E
    if( ((myRank + 1) % root) != 0){
        E = myRank + 1;
    }
 
    //W
    if( (myRank % root) != 0){
        W = myRank - 1;
    }
 
    //NE
    if ( (myRank - root >= 0) && (((myRank - root + 1) % root) != 0)){
        NE = myRank - root + 1;
    }
 
    //NW
    if ( (myRank - root >= 0) && (((myRank - root) % root) != 0)){
        NW = myRank - root - 1;
    }
 
    //SE
    if ( (myRank + root < processes) && (((myRank + root + 1) % root) != 0)){
        SE = myRank + root + 1;
    }
 
    //SW
    if ( (myRank + root < processes) && (((myRank + root) % root) != 0)){
        SW = myRank + root - 1;
    }
    // Create a array for the image section
    unsigned char * grid, * grid2;
    if (color == GRAY){
        grid = calloc( (height/root + 2)*(width/root + 2), sizeof(unsigned char) );   //initialize halo to 0
        if (grid == NULL){
            perror("Calloc");
            return 0;
        }
    }
    else if ( color == COLOR ) {
        grid = calloc( 3*(height/root + 2)*(width/root + 2), sizeof(unsigned char) );   //initialize halo to 0
        if (grid == NULL){
            perror("Calloc");
            return 0;
        }
    }
 
    if (color == GRAY){
        grid2 = calloc( (height/root + 2)*(width/root + 2), sizeof(unsigned char) );   //initialize halo to 0
        if (grid2 == NULL){
            perror("Calloc");
            return 0;
        }
    }
    else if ( color == COLOR ) {
        grid2 = calloc( 3*(height/root + 2)*(width/root + 2), sizeof(unsigned char) );   //initialize halo to 0
        if (grid2 == NULL){
            perror("Calloc");
            return 0;
        }
    }
 
    //open file
    MPI_File fp;
    MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
    int skip = width;
    int gotoMyStart = (myRank/root)*width*(height/root) + (myRank%root)*(width/root);
    unsigned char * readp;
    if (color == GRAY){
        readp = &grid[(width/root + 2 + 1)];
        for (int i=0; i < height/root; i++){
            MPI_File_read_at(fp, gotoMyStart + i*skip, readp, width/root, MPI_BYTE, MPI_STATUS_IGNORE);
            readp = &readp[width/root + 2];
        }
    }
    else if ( color == COLOR ) {
        readp = &grid[3*(width/root + 2 + 1)];
        for (int i=0; i < height/root; i++){
            MPI_File_read_at(fp, 3*gotoMyStart + 3*i*skip, readp, 3*(width/root), MPI_BYTE, MPI_STATUS_IGNORE);
            readp = &readp[3*(width/root + 2)];
        }
    }
    MPI_File_close(&fp);
   
    // MAIN LOOP
    unsigned char * swapGrid = NULL;
    int x = width/root;
    int y = height/root;
    MPI_Barrier(MPI_COMM_WORLD);
    double timer = MPI_Wtime();
   
    for ( int rep = 0; rep < repeats; rep++ ) {
 
        if ( color == GRAY ) {
            // North Data
            MPI_Isend( &grid[ (x+2) + 1 ] , 1 , grayRow, N, 0, MPI_COMM_WORLD, &sendN);
            MPI_Irecv( &grid[ 1 ], 1 , grayRow, N, 0, MPI_COMM_WORLD, &recN);
   
            // North-East Data
            MPI_Isend( &grid[(x+2) + x], 1, grayCorner, NE, 0, MPI_COMM_WORLD, &sendNE);
            MPI_Irecv( &grid[x+1], 1, grayCorner, NE, 0, MPI_COMM_WORLD, &recNE);
   
            // East Data
            MPI_Isend( &grid[ (x+2) + x ] , 1 , grayCol, E, 0, MPI_COMM_WORLD, &sendE);
            MPI_Irecv( &grid[ (x+2) + (x+1) ] , 1 , grayCol, E, 0, MPI_COMM_WORLD, &recE);
   
            // South-East Data
            MPI_Isend( &grid[(x+2)*y + x], 1, grayCorner, SE, 0, MPI_COMM_WORLD, &sendSE);
            MPI_Irecv( &grid[ (x+2)*(y+1) + x+1], 1, grayCorner, SE, 0, MPI_COMM_WORLD, &recSE);
   
            // South Data
            MPI_Isend( &grid[ (x+2)*(y) + 1 ] , 1 , grayRow, S, 0, MPI_COMM_WORLD, &sendS);
            MPI_Irecv( &grid[ (x+2)*(y+1) + 1 ] , 1 , grayRow, S, 0, MPI_COMM_WORLD, &recS);
   
            // South-West Data
            MPI_Isend( &grid[ (x+2)*y + 1 ], 1, grayCorner, SW, 0, MPI_COMM_WORLD, &sendSW);
            MPI_Irecv( &grid[ (x+2)*(y+1) ], 1, grayCorner, SW, 0, MPI_COMM_WORLD, &recSW);
   
            // West Data
            MPI_Isend( &grid[ (x+2) + 1 ] , 1 , grayCol, W, 0, MPI_COMM_WORLD, &sendW);
            MPI_Irecv( &grid[ (x+2) ] , 1 , grayCol, W, 0, MPI_COMM_WORLD, &recW);
   
            // North-West Data
            MPI_Isend( &grid[ (x+2) + 1 ], 1, grayCorner, NW, 0, MPI_COMM_WORLD, &sendNW);
            MPI_Irecv( &grid[ 0 ], 1, grayCorner, NW, 0, MPI_COMM_WORLD, &recNW);
        }
        else if ( color == COLOR ) {
            // North Data
            MPI_Isend( &grid[ 3*(x+2) + 3 ] , 1 , colorRow, N, 0, MPI_COMM_WORLD, &sendN);
            MPI_Irecv( &grid[ 3 ], 1 , colorRow, N, 0, MPI_COMM_WORLD, &recN);
   
            // North-East Data
            MPI_Isend( &grid[ 3*(x+2) + 3*x], 1, colorCorner, NE, 0, MPI_COMM_WORLD, &sendNE);
            MPI_Irecv( &grid[ 3*(x+1) ], 1, colorCorner, NE, 0, MPI_COMM_WORLD, &recNE);
   
            // East Data
            MPI_Isend( &grid[ 3*(x+2) + 3*x ] , 1 , colorCol, E, 0, MPI_COMM_WORLD, &sendE);
            MPI_Irecv( &grid[ 3*(x+2) + 3*(x+1) ] , 1 , colorCol, E, 0, MPI_COMM_WORLD, &recE);
   
            // South-East Data
            MPI_Isend( &grid[ 3*(x+2)*y + 3*x ], 1, colorCorner, SE, 0, MPI_COMM_WORLD, &sendSE);
            MPI_Irecv( &grid[ 3*(x+2)*(y+1) + 3*(x+1) ], 1, colorCorner, SE, 0, MPI_COMM_WORLD, &recSE);
   
            // South Data
            MPI_Isend( &grid[ 3*(x+2)*(y) + 3 ] , 1 , colorRow, S, 0, MPI_COMM_WORLD, &sendS);
            MPI_Irecv( &grid[ 3*(x+2)*(y+1) + 3 ] , 1 , colorRow, S, 0, MPI_COMM_WORLD, &recS);
   
            // South-West Data
            MPI_Isend( &grid[ 3*(x+2)*y + 3 ], 1, colorCorner, SW, 0, MPI_COMM_WORLD, &sendSW);
            MPI_Irecv( &grid[ 3*(x+2)*(y+1) ], 1, colorCorner, SW, 0, MPI_COMM_WORLD, &recSW);
   
            // West Data
            MPI_Isend( &grid[ 3*(x+2) + 3 ] , 1 , colorCol, W, 0, MPI_COMM_WORLD, &sendW);
            MPI_Irecv( &grid[ 3*(x+2) ] , 1 , colorCol, W, 0, MPI_COMM_WORLD, &recW);
   
            // North-West Data
            MPI_Isend( &grid[ 3*(x+2) + 3 ], 1, colorCorner, NW, 0, MPI_COMM_WORLD, &sendNW);
            MPI_Irecv( &grid[ 0 ], 1, colorCorner, NW, 0, MPI_COMM_WORLD, &recNW);
        }
 
       
 
 
        // Convolute the inner image rectangle that isn't related
        // to the other processes data
        convolution(grid,grid2,filter,x+2,2,x-1,2,y-1,color);    
 
 
        // Convolute the rest image using the halo data
        MPI_Wait(&recN, MPI_STATUS_IGNORE);
        MPI_Wait(&recNE, MPI_STATUS_IGNORE);
        MPI_Wait(&recE, MPI_STATUS_IGNORE);
        MPI_Wait(&recSE, MPI_STATUS_IGNORE);
        MPI_Wait(&recS, MPI_STATUS_IGNORE);
        MPI_Wait(&recSW, MPI_STATUS_IGNORE);
        MPI_Wait(&recW, MPI_STATUS_IGNORE);
        MPI_Wait(&recNW, MPI_STATUS_IGNORE);
        // North Data
        convolution(grid,grid2,filter,x+2,2,x-1,1,1,color);
        // North-East Data
        convolution(grid,grid2,filter,x+2,x,x,1,1,color);
        // East Data
        convolution(grid,grid2,filter,x+2,x,x,2,y-1,color);
        // South-East Data
        convolution(grid,grid2,filter,x+2,x,x,y,y,color);
        // South Data
        convolution(grid,grid2,filter,x+2,2,x-1,y,y,color);
        // South-West Data
        convolution(grid,grid2,filter,x+2,1,1,y,y,color);
        // West Data
        convolution(grid,grid2,filter,x+2,1,1,2,y-1,color);
        // North-West Data
        convolution(grid,grid2,filter,x+2,1,1,1,1,color);
 
 
        // Make sure that all sends have finished
        MPI_Wait(&sendN, MPI_STATUS_IGNORE);
        MPI_Wait(&sendNE, MPI_STATUS_IGNORE);
        MPI_Wait(&sendE, MPI_STATUS_IGNORE);
        MPI_Wait(&sendSE, MPI_STATUS_IGNORE);
        MPI_Wait(&sendS, MPI_STATUS_IGNORE);
        MPI_Wait(&sendSW, MPI_STATUS_IGNORE);
        MPI_Wait(&sendW, MPI_STATUS_IGNORE);
        MPI_Wait(&sendNW, MPI_STATUS_IGNORE);
 
        // Output becomes the new input
        swapGrid = grid;
        grid = grid2;
        grid2 = swapGrid;
 
    }
   
 
    timer = MPI_Wtime() - timer;
    double timer2;
    if (myRank != 0){
        MPI_Send(&timer, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else{
        for (int i = 1 ; i < processes ; i++) {
            MPI_Recv(&timer2, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (timer2 > timer)
                timer = timer2;
        }
    }
 
 
    // write to output file
    char *outputFile = malloc((strlen(fileName) + 9) * sizeof(char));
    strcpy(outputFile, "filtered");
    strcat(outputFile, fileName);
    MPI_File_open(MPI_COMM_WORLD, outputFile, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
    if (color == GRAY) {
        readp = &grid[(width/root + 2 + 1)];
        for (int i=0; i < height/root; i++){
            MPI_File_write_at(fp, gotoMyStart + i*skip, readp, width/root, MPI_BYTE, MPI_STATUS_IGNORE);
            readp = &readp[width/root + 2];
        }
    }
    else if ( color == COLOR ) {
            readp = &grid[3*(width/root + 2 + 1)];
        for (int i=0; i < height/root; i++){
            MPI_File_write_at(fp, 3*gotoMyStart + 3*i*skip, readp, 3*(width/root), MPI_BYTE, MPI_STATUS_IGNORE);
            readp = &readp[3*(width/root + 2)];
        }
    }
    MPI_File_close(&fp);
 
    //free
    free(fileName);
    free(outputFile);
    free(grid);
    free(grid2);
   
    MPI_Type_free(&grayRow);
    MPI_Type_free(&colorRow);
    MPI_Type_free(&grayCol);
    MPI_Type_free(&colorCol);
    MPI_Type_free(&grayCorner);
    MPI_Type_free(&colorCorner);
 
    if (myRank == 0){
        printf("Max timer is %f\n", timer);
    }
   
    // Finalize MPI
    MPI_Finalize();
 
}
