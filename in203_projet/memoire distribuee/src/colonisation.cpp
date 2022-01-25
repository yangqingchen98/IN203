#include <cstdlib>
#include <string>
#include <iostream>
#include <SDL.h>
#include <SDL_image.h>
#include <fstream>
#include <ctime>
#include <iomanip>      // std::setw
#include <chrono>
#include <mpi.h>

#include "parametres.hpp"
#include "galaxie.hpp"
 
int main(int argc, char ** argv)
{
    
    char commentaire[4096];
    int width, height, nbp,rank, height_pas;
    
    
    parametres param;

    MPI_Init(&argc, &argv);
    MPI_Comm globComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
    MPI_Comm_size(globComm, &nbp);
    MPI_Comm_rank(globComm, &rank);
    MPI_Status status;
    
    std::ifstream fich("parametre.txt");
    fich >> width;
    fich.getline(commentaire, 4096);
    fich >> height;
    fich.getline(commentaire, 4096);
    fich >> param.apparition_civ;
    fich.getline(commentaire, 4096);
    fich >> param.disparition;
    fich.getline(commentaire, 4096);
    fich >> param.expansion;
    fich.getline(commentaire, 4096);
    fich >> param.inhabitable;
    fich.getline(commentaire, 4096);
    fich.close();

    height_pas = int(height/(nbp-1));
    
    if (rank == 0){
        
        SDL_Event event;
        SDL_Window   * window;
        
        std::cout << "Resume des parametres (proba par pas de temps): " << std::endl;
        std::cout << "\t Chance apparition civilisation techno : " << param.apparition_civ << std::endl;
        std::cout << "\t Chance disparition civilisation techno: " << param.disparition << std::endl;
        std::cout << "\t Chance expansion : " << param.expansion << std::endl;
        std::cout << "\t Chance inhabitable : " << param.inhabitable << std::endl;
        std::cout << "Proba minimale prise en compte : " << 1./RAND_MAX << std::endl;
        std::srand(std::time(nullptr));
		
        SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO);
        window = SDL_CreateWindow("Galaxie", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                      width, height, SDL_WINDOW_SHOWN);
        
        galaxie g(width, height, param.apparition_civ);
        galaxie_renderer gr(window);
        
        int deltaT = (20*52840)/width;
        std::cout << "Pas de temps : " << deltaT << " années" << std::endl;

        std::cout << std::endl;
		
        gr.render(g);
        unsigned long long temps = 0;
            
        std::vector<char> S_Buffer(width*height_pas);
        std::vector<char> R_Buffer(width*height_pas);
        std::vector<char> S_Buffer_Edge(width);
            
        std::chrono::time_point<std::chrono::system_clock> start, end1, end2;
        for(int count = 1; count < nbp; ++count){
            for(int i = 0; i < height_pas; ++i){
                for(int j = 0; j < width; ++j){
                    S_Buffer[i*width+j] = *(g.data()+(count-1)*width*height_pas+i*width+j);
                }
            }
            MPI_Send(S_Buffer.data(), width*height_pas, MPI_CHAR, count, 0, globComm);
        }
        
        while(1){
        	start = std::chrono::system_clock::now();
        	for(int count = 1; count < nbp-1; ++count){
                for(int j = 0; j < width; ++j){
                    S_Buffer_Edge[j] = *(g.data()+count*width*height_pas+j);
                }
                MPI_Send(S_Buffer_Edge.data(), width, MPI_CHAR, count, 0, globComm);
            }
            
            for(int count = 2; count < nbp; ++count){
                for(int j = 0; j < width; ++j){
                    S_Buffer_Edge[j] = *(g.data()+(count-1)*width*height_pas-width+j);   
                }
                MPI_Send(S_Buffer_Edge.data(), width, MPI_CHAR, count, 0, globComm);
            }
            
            for(int count = 1; count < nbp; ++count){
                MPI_Recv(R_Buffer.data(), width*height_pas, MPI_CHAR, count, 0, globComm, &status);
                g.SetValue(R_Buffer, (count-1)*width*height_pas, height_pas);
            }
            
            end1 = std::chrono::system_clock::now();
            gr.render(g);
            
            end2 = std::chrono::system_clock::now();

                std::chrono::duration<double> elaps1 = end1 - start;
                std::chrono::duration<double> elaps2 = end2 - end1;
                
                temps += deltaT;
                std::cout << "Temps passe : "
                          << std::setw(10) << temps << " années"
                          << std::fixed << std::setprecision(3)
                          << "  " << "|  CPU(ms) : calcul " << elaps1.count()*1000
                          << "  " << "affichage " << elaps2.count()*1000
                          << "\r" << std::flush;
                //_sleep(1000);
             
            if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                std::cout << std::endl << "The end" << std::endl;
                SDL_DestroyWindow(window);
    			SDL_Quit();
    			MPI_Finalize();
                break;
            }
        }
    }else if(rank == 1){
        std::vector<char> S_Buffer(width*height_pas);
        std::vector<char> R_Buffer(width*height_pas);
        std::vector<char> R_Buffer_Edge(width);
        std::vector<char> g(width*height_pas);
        std::vector<char> g_next(width*height_pas);
        MPI_Recv(R_Buffer.data(), width*height_pas, MPI_CHAR, 0, 0, globComm, &status);
		g.assign(R_Buffer.begin(), R_Buffer.end());
        while(1){
        	MPI_Recv(R_Buffer_Edge.data(), width, MPI_CHAR, 0, 0, globComm, &status);
            g.insert(g.end(), R_Buffer_Edge.begin(), R_Buffer_Edge.end());

            mise_a_jour(param, width, height_pas+1, g.data(), g_next.data());

            S_Buffer.assign(g_next.begin(), g_next.end()-width);
            g.assign(S_Buffer.begin(), S_Buffer.end());
            MPI_Send (S_Buffer.data(), width*height_pas, MPI_CHAR, 0, 0, globComm);
        }
    }else if(rank == nbp-1){
    	std::vector<char> S_Buffer(width*height_pas);
        std::vector<char> R_Buffer(width*height_pas);
        std::vector<char> R_Buffer_Edge(width);
        std::vector<char> g(width*height_pas);
        std::vector<char> g_next(width*height_pas);
        MPI_Recv(R_Buffer.data(), width*height_pas, MPI_CHAR, 0, 0, globComm, &status);
        g.assign(R_Buffer.begin(), R_Buffer.end());
        while(1){
        	MPI_Recv(R_Buffer_Edge.data(), width, MPI_CHAR, 0, 0, globComm, &status);
            g.insert(g.begin(), R_Buffer_Edge.begin(), R_Buffer_Edge.end());
            
            mise_a_jour(param, width, height_pas+1, g.data(), g_next.data());

            S_Buffer.assign(g_next.begin()+width, g_next.end());
            g.assign(S_Buffer.begin(), S_Buffer.end());
            MPI_Send (S_Buffer.data(), width*height_pas, MPI_CHAR, 0, 0, globComm);
        }
    }else{
    	std::vector<char> S_Buffer(width*height_pas);
        std::vector<char> R_Buffer(width*height_pas);
        std::vector<char> R_Buffer_Edge(width);
        std::vector<char> g(width*height_pas);
        std::vector<char> g_next(width*height_pas);
        MPI_Recv(R_Buffer.data(), width*height_pas, MPI_CHAR, 0, 0, globComm, &status);
        g.assign(R_Buffer.begin(), R_Buffer.end());
        while(1){
        	MPI_Recv(R_Buffer_Edge.data(), width, MPI_CHAR, 0, 0, globComm, &status);
            g.insert(g.end(), R_Buffer_Edge.begin(), R_Buffer_Edge.end());
            MPI_Recv(R_Buffer_Edge.data(), width, MPI_CHAR, 0, 0, globComm, &status);
            g.insert(g.begin(), R_Buffer_Edge.begin(), R_Buffer_Edge.end());

            mise_a_jour(param, width, height_pas+2, g.data(), g_next.data());

            S_Buffer.assign(g_next.begin()+width, g_next.end()-width);
            g.assign(S_Buffer.begin(), S_Buffer.end());
            MPI_Send (S_Buffer.data(), width*height_pas, MPI_CHAR, 0, 0, globComm);
        }
    }
    
    std::cout << rank << " MPI_Finalize "<< std::endl;
    return EXIT_SUCCESS;
}
