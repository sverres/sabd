//
//       name       : sabd.cpp
//       version    : 1.0
//       date       : 99.05.07
//       made by    : sverre stikbakke
//       description: Simulated annealing algorithm used
//                    on biodiversity data.
//
//       syntax     :  sabd      ck_0   Lk_0   red_rate   outlevel
//
//       reference  : Aarts and Korst, 1989.
//                    Simulated Annealing and Boltzmann Machines.
//


#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SPECIES 133          //number of species
#define UNITS 155            //number of UNITS
#define MAX_SELECTED 10      //number of UNITS to be selected

int data[UNITS][SPECIES] =   //matrix of species composition
#include "species.txt"

int solution[MAX_SELECTED];               //vector of UNITS currently in the solution
int candidate_solution[MAX_SELECTED];     //vector of UNITS to be evaluated as the next solution

float ck;               //control parameter (temperature)
int Lk;                 //step-length (iterations on each temperature)

float ck_0;             // initial temperature
int   Lk_0;             // initial step-length
float red_rate;         // temperature reduction rate
int   outlevel;         // output level control

int k;                  //step counter
int l;                  //iteration counter at each step
int m;                  //iteration counter
int no_improvement;     //counter for non-improving steps
int accepted;           //counter for accepted solutions on each step
int Z;                  //objective function value
int new_Z;              //objective function value of new solution
int max_Z;              //max value on each step

// function prototypes

void print_usage();
void init(int * _solution);
int  stopcriterion();
void generate(int * _solution, int * _candidate_solution);
int  eval_obj_f(int * _canditate_solution);
void print_solution( int * _solution, int _Z );
int  accept( int _new_Z, int _Z, int * _candidate_solution, int * _solution);
void insert( int * _candidate_solution, int * _solution );
void calculate_length();
void calculate_control(float _red_rate);
void print_result();
int  random(int _seed, int _mod);


// main function

int main( int argc, char *argv[])
{

    // command line argument handling

    if (argc < 4)
    {
        print_usage();

        return -1;
    }

    sscanf(argv[ 1 ],"%g", &ck_0);
    Lk_0 = atoi(argv[ 2 ]);
    sscanf(argv[ 3 ],"%g", &red_rate);
    outlevel = atoi(argv[ 4 ]);

    // initialize parameters and make first solution

    ck = ck_0;
    Lk = Lk_0;

    init(solution);

    Z = eval_obj_f(solution);

    max_Z = Z;

    if ( outlevel > 5 ) print_solution( solution, Z );

    // main loop

    while (!stopcriterion())
    {
        // loop on each temperature

        max_Z = 0;

        for ( int l = 0; l < Lk; l++ )
        {
            generate(solution, candidate_solution);

            new_Z = eval_obj_f(candidate_solution);

            //save max Z value on each step
            if (new_Z > max_Z) {max_Z = new_Z;}

            if ( outlevel > 6 ) print_solution( candidate_solution, new_Z );
            if (( outlevel > 5 ) && (l > (Lk - 5 ))) print_solution( candidate_solution, new_Z );

            if (accept(new_Z, Z, candidate_solution, solution))
            {
                if ( outlevel == 5 ) print_solution( candidate_solution, new_Z);
                if ( outlevel == 7 ) cout << " - - - solution accepted\n";

                insert( candidate_solution, solution );

                Z = new_Z;

                accepted++;
            }

            m++;        //increment iteration counter
        }

        k++;        //increment step counter


        // print results from each step

        if ( outlevel == 1 ) cout << ".";
        if ( outlevel > 4 ) cout << "\n";
        if ( outlevel > 1 )
        {
             cout << "Max_Z: " << max_Z << "   step: " << k+1 << "   ck: " << ck
                  << "\    accepted/trials: "<< accepted << "/" << Lk << "\n";
        }
        if ( outlevel > 4 ) cout << "\n";

        // adjust parameters

        if (accepted == 0)
        {
            no_improvement++;
        }
        else
        {
            no_improvement = 0;
        }

        accepted = 0;

        calculate_length();  //dummy

        calculate_control( red_rate );
    }

    print_result();

    return 0;
}


// init - set step length, control parameter
// and a random initial soluton

void init(int * _solution )
{

    int seed =3412;

    for (int i = 0; i < MAX_SELECTED; i++)
    {
        _solution[ i ] = random( seed, UNITS );

        seed = seed + 173;
    }
}


// generate - make a new solution by changing one unit randomly

void generate(int * _solution, int * _candidate_solution)
{
    int seed = 0;

    for (int i = 0; i < MAX_SELECTED; i++)
    {
        _candidate_solution[ i ] = _solution[ i ];
    }

    seed = m * 50 + 7369;  // arbritary chosen value
    int position_to_change = random( seed, MAX_SELECTED );

    seed = m * 10 + 3521;  // arbritary chosen value
    int new_unit_selected  = random( seed, UNITS);

    _candidate_solution[ position_to_change ] = new_unit_selected;

/*

    cout << "k :"  << k;
    cout << "\tchange_index:  " << change_index;
    cout << "\tnew_candidate: " << new_candidate << "\n";

    cout << "Old: ";

    for ( i = 0; i < MAX_SELECTED; i++)
    {
        cout << _solution[ i ] << "  ";
    }

    cout << "\nNew: ";

    for ( i = 0; i < MAX_SELECTED; i++)
    {
        cout << _candidate_solution[ i ] << "  ";
    }

    cout << "\n";
*/

}


// eval_obj_f - counts species in solution

int eval_obj_f(int * _solution)
{
    int units_with_specie = 0;
    int species_in_solution = 0;

    // loop over all species

    for (int i = 0; i < SPECIES; i++)
    {

        // loop over all units in _solution
        // and count actual specie

        units_with_specie = 0;

        for (int j = 0; j < MAX_SELECTED; j++)
        {
          units_with_specie = units_with_specie + data[_solution[j]][i];
        }

        if (units_with_specie > 0)
        {
            species_in_solution++;
        }

    }

    return species_in_solution;
}


// print_solution

void print_solution( int * _solution, int _Z )
{

    cout << "Z = " << _Z << "    ";
    for ( int i = 0; i < MAX_SELECTED; i++)
    {
        //prints unit number starting with index 1
        cout << _solution[ i ] + 1  << "  ";
    }

    cout << "    m: " << m << "\n";
}



// accept - decides whether the candidate solution will
//          replace current solution

int accept( int _new_Z, int _Z, int * _candidate_solution, int * _solution)
{
    int yes_no  = 0;
    int diff    = 0;

    float eval = 0;
    float p    = 0;

    // test if the _cadidate_solution differs from _solution

    for (int i = 0; i < MAX_SELECTED; i++)
    {
        if (_solution[ i ] != _candidate_solution[ i ])
        {
            diff++;
        }
    }

    if (diff == 0)
    {
        yes_no = 0;
        return yes_no;
    }

    // evaluate new_Z against Z under
    // influence of temperature ck

    if (new_Z >= Z)
    {
        yes_no = 1;
        eval = 1;
        p = 1;
    }
    else
    {

        float A    = (_new_Z - _Z) / ck;
        eval = exp( A );

        int   seed = m * 13;

        //generate random number in [0, 1]

        p    = random( seed, 1000) * 0.001;

        // accept a non-improving solution?

        if (eval > p ) { yes_no = 1; }
    }

/*
    if (_l < 5 || _l > (Lk-5) )
    {
        cout << "new_Z: " << _new_Z << "\taccept: " << yes_no << "\tp: "
        << p << "\teval: " << eval << "\tck: " << ck << "\tk: " << k << "\n";
    }
*/
    return yes_no;
}

// insert new solution

void insert( int * _candidate_solution, int * _solution )
{

    for (int i = 0; i < MAX_SELECTED; i++)
    {
        _solution[ i ] = _candidate_solution[ i ];
    }

}


// calculate step-length (dummy)

void calculate_length(){}

// calculate_control: lower temperature

void calculate_control(float _red_rate)
{
    ck = ck - ck * _red_rate;
}

// stopcriterion based on non-improving steps

int stopcriterion()
{
    if (no_improvement > 5 || ck < 0.00001 )
    {
        return 1;
    }
    else return 0;
}

// random number function

int random (int _seed, int _max_val)
{
    srand(_seed);

    return (rand() % _max_val);
}

void print_usage()
{
    cout << "\nUsage: \n\n";

    cout << "sabd      ck_0   Lk_0   red_rate   outlevel\n\n";

    cout << "where\n\n";

    cout << "ck_0     = initial temperature             (e.g. 35)  \n";
    cout << "Lk_0     = iterations on each step         (e.g. 50)  \n";
    cout << "red_rate = temperature reduction rate     (e.g. 0.2)  \n";
    cout << "outlevel = output level control        (0,1,2,5,6,7)  \n";
    cout << "               0 : only summary                                    \n";
    cout << "               1 : one dot on each step                            \n";
    cout << "               2 : summary data from each step                     \n";
    cout << "               5 : prints every accepted solution                  \n";
    cout << "               6 : prints 5 last generated solutions on each step  \n";
    cout << "               7 : prints all generated solutions                  \n";
    cout << "\n";
}



void print_result()
{
    cout << "\n";

    cout << "############## Simulated annealing results summary ##################\n\n";

    cout << "Parameters: \n\n";

    cout << "initial temperature,     ck_0: " << ck_0           << "\n";
    cout << "iterations on each step, Lk_0: " << Lk_0           << "\n";
    cout << "temperature reduction rate   : " << red_rate       << "\n";
    cout << "output level control         : " << outlevel       << "\n";
    cout << "\n";

    cout << "Results: \n\n";

    cout << "Steps (k):               " << k  << "\n";
    cout << "Iterations (m):          " << m  << "\n";
    cout << "Final state temperature: " << ck << "\n";
    cout << "\n";

    cout << "Last accepted solution: \n\n";

    print_solution( solution, Z );

    cout << "\n#####################################################################\n";

}
