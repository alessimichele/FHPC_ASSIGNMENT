#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>


void ordered_update_finite(unsigned char* grid, int k, int n_steps, int s, int ){
    /*
    evolve the current state of the game of life grid for n_steps using the ordered 
    update algorithm. The grid has size k x k.
    Parameters
    ----------
    grid: pointer to unsigned char, the grid, which is a 1D array, where 
        grid[i*k+j] is the state of the cell in row i and column j
    k: int, size of the grid
    n_steps: int, number of steps to evolve the grid
    s: int, every how many steps a dump of the system is saved on a file
        (0 meaning only at the end)
    */
    int n_threads = 1   
    unsigned char* flag = (unsigned char*)malloc(k*k*sizeof(unsigned char));
    #pragma omp parallel
    {   
        my_id = omp_get_thread_num();
        #pragma omp master{
            int n_threads = omp_get_num_threads();
            if (n_threads > k)
            {
                omp_set_num_threads(k);
                n_threads = k;
            }
        }
        //now we want to update the matrix so that each thread updates a row, waiting for the previous
        //row to be two steps ahead
        //lines_per_thread is going to be the floor of k/n_threads
        int rows_for_me = k/n_threads; //floor of k/n_threads
        if (my_id < k%n_threads)
        {
            rows_for_me++; //total number of rows to be updated by each thread
        }
        for (int step = 0; step < n_steps; step++)
        {   for (int line = 0; line < rows_for_me, line++)
            {   
                int my_current_row = my_id + n_threads * line; //my_id is the rest of the division of k by n_threads
                if (my_current_row==0)//first row, different update since grid is not infinite
                {    //first cell
                    int n_neigh_255 = grid[1] + grid[k] + grid[k+1];
                    grid[0] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                    flag[0] = 1;
                    //middle cells
                    for (int j=1; j<k-1; j++)
                    {             
                        n_neigh_255 = grid[j-1] + 
                                      grid[j+1] + 
                                      grid[k+j-1] + 
                                      grid[k+j] + 
                                      grid[k+j+1];
                        grid[j] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                        flag[j] = 1;
                    }
                    //no need for while for the last cell of the row
                    n_neigh_255 = grid[k-2] + grid[2*k-2] + grid[2*k-1];
                    grid[k-1] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                    flag[k-1] = 1;
                }else if (my_current_row==k-1)//last row, different update since grid is not infinite
                {//first cell
                    while (flag[(k-2)*k+1] !=1){}  
                    n_neigh_255 = grid[(k-2)*k] + grid[(k-2)*k+1] + grid[(k-1)*k+1];
                    grid[(k-1)*k] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                    flag[(k-1)*k] = 1;
                    //middle cells
                    for (int j=1; j<k-1; j++)
                    {   
                        while (flag[(k-2)*k+j+1] !=1){}   //wait for the previous thread to update two cells in the row above           
                        n_neigh_255 = grid[(k-1)*k+j-1] + 
                                      grid[(k-1)*k+j+1] + 
                                      grid[(k-2)*k+j-1] + 
                                      grid[(k-2)*k+j] + 
                                      grid[(k-2)*k+j+1];
                        grid[(k-1)*k+j] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                        flag[(k-1)*k+j] = 1;
                    }
                    //no need for while for the last cell of the row
                    n_neigh_255 = grid[(k-2)*k+k-2] + grid[(k-2)*k+k-1] + grid[(k-1)*k+k-2];
                    grid[(k-1)*k+k-1] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                    //flag[(k-1)*k+k-1] = 1; no need, it's the last one
                }else
                {    
                    while (flag[(my_current_row-1)*k+1] !=1){}   //wait for the previous thread to update two cells in the row above
                    n_neigh_255 = grid[(my_current_row-1)*k] + 
                                  grid[(my_current_row-1)*k+1] +
                                  grid[my_current_row*k+1] + 
                                  grid[(my_current_row+1)*k] +
                                  grid[(my_current_row+1)*k+1];
                    grid[my_current_row*k] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                    flag[my_current_row*k] = 1;
                    for (int j=1; j<k-1; j++)
                    {   
                        while (flag[(my_current_row-1)*k+j+1] !=1){}   //wait for the previous thread to update two cells in the row above           
                        n_neigh_255 = grid[my_current_row*k+j-1] + 
                                      grid[my_current_row*k+j+1] + 
                                      grid[(my_current_row-1)*k+j-1] + 
                                      grid[(my_current_row-1)*k+j] + 
                                      grid[(my_current_row-1)*k+j+1]+
                                      grid[(my_current_row+1)*k+j-1] +
                                      grid[(my_current_row+1)*k+j] +
                                      grid[(my_current_row+1)*k+j+1];
                        grid[i*k+j] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                        flag[my_current_row*k+j] = 1;
                    }
                    //no need for while for the last cell of the row
                    n_neigh_255 = grid[(my_current_row-1)*k+k-2] + 
                                  grid[(my_current_row-1)*k+k-1] + 
                                  grid[my_current_row*k+k-2] + 
                                  grid[(my_current_row+1)*k+k-2] +
                                  grid[(my_current_row+1)*k+k-1];
                    grid[my_current_row*k+k-1] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
                    flag[my_current_row*k+k-1] = 1;
                }
            }
            #pragma omp barrier
            #pragma omp master{
                if (step+1%s==0)
                {
                    char path[45] = "images/evolve_ordered/";
                    char name[20];
                    snprintf(name, 20, "snapshot_%05d.pgm", current_step+1);
                    strcat(path, name);
                    //DA RIVEDERE QUANDO AVREMO SCRITTO LA FUNZIONE write_pgm_parallel
                    write_pgm_parallel(grid, 255, k, k, path, 0, 1, k);
                }
            }
        }
    }
    free(flag);
}

/*generica idea di documentazione: siccome nella grid finita non è possibile
fare i trucchetti col modulo, siamo costretti a fare un update ordinato, cioè
caso per caso (prima riga, ultima riga, righe intermedie). In questo modo tutti
i bordi sono aggiornati per conto loro.
Il while fa un po' schifo, ma non ho trovato un modo migliore interno a openmp
per farlo funzionare. In pratica, ogni thread aspetta a vuoto che il thread precedente
abbia aggiornato due celle nella riga sopra prima di aggiornare la sua riga. 
Il modo in cui controllo che il thread precedente abbia aggiornato due celle è
tramite un array flag della stessa dimensione della grid, le cui celle sono 
settate a 1 quando il thread ha aggiornato la cella corrispondente.
*/