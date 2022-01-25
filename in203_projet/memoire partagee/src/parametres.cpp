#include <cstdlib>
#include <cassert>
#include <ctime>
#include <iostream>
#include <chrono>
#include <stdlib.h>
#include <random>
#include "galaxie.hpp"
#include "parametres.hpp"
#include "omp.h"


expansion calcul_expansion(const parametres& c, std::uniform_real_distribution<double>& urd,
                           std::minstd_rand0& r0)
{
    double val = urd(r0);
    if (val < 0.01*c.expansion)     // parmi c.expansion, on a 1% de chance d'expansion isotrope...
        return expansion_isotrope;
    if (val < c.expansion)          // ... et 99% de chance d'expansion dans 1 seule direction
        return expansion_unique;
    return pas_d_expansion;
}
//_ ______________________________________________________________________________________________ _
bool calcul_depeuplement(const parametres& c, std::uniform_real_distribution<double>& urd,
                         std::minstd_rand0& r0)
{
    double val = urd(r0);
    if (val < c.disparition)
        return true;
    return false;   
}
//_ ______________________________________________________________________________________________ _
bool calcul_inhabitable(const parametres& c, std::uniform_real_distribution<double>& urd,
                        std::minstd_rand0& r0)
{
    double val = urd(r0);
    if (val < c.inhabitable)
        return true;
    return false;
}
//_ ______________________________________________________________________________________________ _
bool apparition_technologie(const parametres& p, std::uniform_real_distribution<double>& urd,
                            std::minstd_rand0& r0)
{
    double val = urd(r0);
    if (val < p.apparition_civ)
        return true;
    return false;
}
//_ ______________________________________________________________________________________________ _
bool a_un_systeme_proche_colonisable(int i, int j, int width, int height, const char* galaxie)
{
    assert(i >= 0);
    assert(j >= 0);
    assert(i < height);
    assert(j < width);

    if ( (i>0) && (galaxie[(i-1)*width+j] == habitable)) return true;
    if ( (i<height-1) && (galaxie[(i+1)*width+j] == habitable)) return true;
    if ( (j>0) && (galaxie[i*width+j-1] == habitable)) return true;
    if ( (j<width-1) && (galaxie[i*width+j+1] == habitable)) return true;

    return false;
}
//_ ______________________________________________________________________________________________ _

void 
mise_a_jour(const parametres& params, int width, int height, const char* galaxie_previous, char* galaxie_next)
{
    int i, j;
    
    memcpy(galaxie_next, galaxie_previous, width*height*sizeof(char));
    
    #pragma omp parallel for private(i,j) shared(galaxie_previous, galaxie_next)
    for ( i = 0; i < height; ++i )
      {
        for ( j = 0; j < width; ++j )
        {
            unsigned graine = std::chrono::system_clock::now().time_since_epoch().count() + 18679*omp_get_thread_num();
            std::minstd_rand0 r0(graine);
            std::uniform_real_distribution<double> urd(0.0,1.0);
            std::uniform_int_distribution<int> uid(0,3);

            if (galaxie_previous[i*width+j] == habitee)
            {
                if ( a_un_systeme_proche_colonisable(i, j, width, height, galaxie_previous) )
                {
                    expansion e = calcul_expansion(params, urd, r0);
                    if (e == expansion_isotrope)
                    {
                      if ( (i > 0) && (galaxie_previous[(i-1)*width+j] != inhabitable) )
                        {
                            galaxie_next[(i-1)*width+j] = habitee;
                        }
                      if ( (i < height-1) && (galaxie_previous[(i+1)*width+j] != inhabitable) )
                        {
                            galaxie_next[(i+1)*width+j] = habitee;
                        }
                      if ( (j > 0) && (galaxie_previous[i*width+j-1] != inhabitable) )
                        {
                            galaxie_next[i*width+j-1] = habitee;
                        }
                      if ( (j < width-1) && (galaxie_previous[i*width+j+1] != inhabitable) )
                        {
                            galaxie_next[i*width+j+1] = habitee;
                        }
                    }
                    else if (e == expansion_unique)
                    {
                        // Calcul de la direction de l'expansion :
                        int ok = 0;
                        do
                        {
                            int dir = uid(r0);
                            if ( (i>0) && (0 == dir) && (galaxie_previous[(i-1)*width+j] != inhabitable) )
                            {
                                galaxie_next[(i-1)*width+j] = habitee;
                                ok = 1;
                            }
                            if ( (i<height-1) && (1 == dir) && (galaxie_previous[(i+1)*width+j] != inhabitable) )
                            {
                                galaxie_next[(i+1)*width+j] = habitee;
                                ok = 1;
                            }
                            if ( (j>0) && (2 == dir) && (galaxie_previous[i*width+j-1] != inhabitable) )
                            {
                                galaxie_next[i*width+j-1] = habitee;
                                ok = 1;
                            }
                            if ( (j<width-1) && (3 == dir) && (galaxie_previous[i*width+j+1] != inhabitable) )
                            {
                                galaxie_next[i*width+j+1] = habitee;
                                ok = 1;
                            }
                        } while (ok == 0);
                    }// End if (e == expansion_unique)
                }// Fin si il y a encore un monde non habite et habitable
                if (calcul_depeuplement(params, urd, r0))
                {
                    galaxie_next[i*width+j] = habitable;
                }
                if (calcul_inhabitable(params, urd, r0))
                {
                    galaxie_next[i*width+j] = inhabitable;
                }
            }  // Fin si habitee
            else if (galaxie_previous[i*width+j] == habitable)
            {
                if (apparition_technologie(params, urd, r0))
                    galaxie_next[i*width+j] = habitee;
            }
            else { // inhabitable
              // nothing to do : le systeme a explose
            }
            // if (galaxie_previous...)
        }// for (j)
      }// for (i)

}
//_ ______________________________________________________________________________________________ _
