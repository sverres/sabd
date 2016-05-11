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

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SPECIES 133     // number of species
#define UNITS 155       // number of UNITS
#define MAX_SELECTED 10 // number of UNITS to be selected

#define IMPROVEMENT_CRITERION 5       // part of stopcriterion
#define TEMPERATURE_CRITERION 0.00001 // part of stopcriterion

int data[UNITS][SPECIES] = // matrix of species composition
#include "species.txt"

int solution[MAX_SELECTED]; // vector of UNITS currently in the solution
int candidate_solution[MAX_SELECTED]; // vector of UNITS to be evaluated as the
                                      // next solution

float ck; // control parameter (temperature)
int Lk;   // step-length (iterations on each temperature)

float ck_0;     // initial temperature
int Lk_0;       // initial step-length
float red_rate; // temperature reduction rate
int outlevel;   // output level control

int k;              // step counter
int l;              // iteration counter at each step
int m;              // iteration counter
int no_improvement; // counter for non-improving steps
int accepted;       // counter for accepted solutions on each step
int Z;              // objective function value
int new_Z;          // objective function value of new solution
int max_Z;          // max value on each step

// function prototypes

void print_usage();
int random(int _max_val);
void init(int *_solution);
bool unique_units(int *_solution);
bool stopcriterion();
void generate(int *_solution, int *_candidate_solution);
int eval_obj_f(int *_canditate_solution);
bool accept(int _new_Z, int _Z, int *_candidate_solution, int *_solution);
void insert(int *_candidate_solution, int *_solution);
void calculate_control(float _red_rate);
void print_solution(int *_solution, int _Z);
void print_result();

using std::cout;

//----------------------------------------------------------------------------//
// main function
//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {

  // command line argument handling

  if (argc < 4) {
    print_usage();
    return -1;
  }

  sscanf(argv[1], "%g", &ck_0);
  Lk_0 = atoi(argv[2]);
  sscanf(argv[3], "%g", &red_rate);
  outlevel = atoi(argv[4]) || 0;

  // initialize parameters and make first solution

  ck = ck_0;
  Lk = Lk_0;

  srand(time(NULL)); // random number generator initialization

  do {
    init(solution);
  } while (!unique_units(solution));

  Z = eval_obj_f(solution);

  max_Z = Z;

  // main loop

  while (!stopcriterion()) {

    // loop on each step (temperature)

    max_Z = 0;

    for (int l = 0; l < Lk; l++) {
      do {
        generate(solution, candidate_solution);
      } while (!unique_units(candidate_solution));

      new_Z = eval_obj_f(candidate_solution);

      // save max Z value on each iteration

      if (new_Z > max_Z) {
        max_Z = new_Z;
      }

      if (accept(new_Z, Z, candidate_solution, solution)) {
        insert(candidate_solution, solution);
        Z = new_Z;
        accepted++;
      }
      m++; // iteration counter
    }
    k++; // step counter

    // update control parameters

    if (accepted == 0) {
      no_improvement++;
    } else {
      no_improvement = 0;
    }

    accepted = 0;
    calculate_control(red_rate);

  } // end main loop

  print_result();
  return 0;
}

//----------------------------------------------------------------------------//
// init - generate a random initial soluton
//----------------------------------------------------------------------------//
void init(int *_solution) {
  for (int i = 0; i < MAX_SELECTED; i++) {
    _solution[i] = random(UNITS);
  }
}

//----------------------------------------------------------------------------//
// test if units in solution are unique
//----------------------------------------------------------------------------//
bool unique_units(int *_solution) {
  for (int i = 0; i < MAX_SELECTED; i++) {
    for (int j = 0; j < MAX_SELECTED; j++) {
      if ((i != j) && (_solution[i] == _solution[j])) {
        return false;
      };
    }
  }
  return true;
}

//----------------------------------------------------------------------------//
// generate - make a new solution by changing one unit randomly
//----------------------------------------------------------------------------//
void generate(int *_solution, int *_candidate_solution) {
  for (int i = 0; i < MAX_SELECTED; i++) {
    _candidate_solution[i] = _solution[i];
  }

  int position_to_change = random(MAX_SELECTED);
  int new_unit_selected = random(UNITS);

  _candidate_solution[position_to_change] = new_unit_selected;
}

//----------------------------------------------------------------------------//
// eval_obj_f - counts species in solution
//----------------------------------------------------------------------------//
int eval_obj_f(int *_solution) {
  int units_with_specie = 0;
  int species_in_solution = 0;

  // loop over all species

  for (int i = 0; i < SPECIES; i++) {

    // loop over all units in _solution
    // and count actual specie

    units_with_specie = 0;
    for (int j = 0; j < MAX_SELECTED; j++) {
      units_with_specie = units_with_specie + data[_solution[j]][i];
    }

    if (units_with_specie > 0) {
      species_in_solution++;
    }
  }
  return species_in_solution;
}

//----------------------------------------------------------------------------//
// print_solution
//----------------------------------------------------------------------------//
void print_solution(int *_solution, int _Z) {
  cout << "Z = " << _Z << "    ";
  for (int i = 0; i < MAX_SELECTED; i++) {

    // prints unit number starting with index 1
    cout << _solution[i] + 1 << "  ";
  }
  cout << "    m: " << m << "\n";
}

//----------------------------------------------------------------------------//
// accept - decides whether the candidate solution will
//          replace current solution
//----------------------------------------------------------------------------//
bool accept(int _new_Z, int _Z, int *_candidate_solution, int *_solution) {

  int diff = 0;
  for (int i = 0; i < MAX_SELECTED; i++) {
    if (_solution[i] != _candidate_solution[i]) {
      diff++;
    }
  }

  if (diff == 0) {
    return false;
  }

  if (new_Z >= Z) {
    return true;
  } else {

    // accept a non-improving solution?

    // evaluate _new_Z against _Z under
    // influence of temperature ck

    float A = (_new_Z - _Z) / ck;
    float eval = exp(A);

    // generate random number in [0, 1]
    float p = random(1000) * 0.001;

    if (eval > p) {
      return true;
    }
  }
  return false;
}

//----------------------------------------------------------------------------//
// insert new solution
//----------------------------------------------------------------------------//
void insert(int *_candidate_solution, int *_solution) {
  for (int i = 0; i < MAX_SELECTED; i++) {
    _solution[i] = _candidate_solution[i];
  }
}

//----------------------------------------------------------------------------//
// calculate_control: lower temperature
//----------------------------------------------------------------------------//
void calculate_control(float _red_rate) { ck = ck - ck * _red_rate; }

//----------------------------------------------------------------------------//
// stopcriterion based on non-improving steps
//----------------------------------------------------------------------------//
bool stopcriterion() {
  if (no_improvement > IMPROVEMENT_CRITERION || ck < TEMPERATURE_CRITERION) {
    return true;
  } else
    return false;
}

//----------------------------------------------------------------------------//
// random number function (seeded in main function)
//----------------------------------------------------------------------------//
int random(int _max_val) { return (rand() % _max_val); }

//----------------------------------------------------------------------------//
// print usage
//----------------------------------------------------------------------------//
void print_usage() {
  cout << "\nUsage: \n\n";

  cout << "sabd      ck_0   Lk_0   red_rate   outlevel\n\n";
  cout << "where\n\n";

  cout << "ck_0     = initial temperature             (e.g. 35)   \n";
  cout << "Lk_0     = iterations on each step         (e.g. 50)   \n";
  cout << "red_rate = temperature reduction rate     (e.g. 0.2)   \n";
  cout << "outlevel = output level control        (0,1,2,5,6,7) \n\n";
  cout << "  0 : only summary                                     \n";
  cout << "  1,2,5,6,7 : not in effect in this version            \n";
}

//----------------------------------------------------------------------------//
// print results
//----------------------------------------------------------------------------//
void print_result() {
  cout << "\n############## Simulated annealing results summary ###############"
          "###\n\n";

  cout << "Parameters: \n\n";

  cout << "initial temperature,     ck_0: " << ck_0 << "\n";
  cout << "iterations on each step, Lk_0: " << Lk_0 << "\n";
  cout << "temperature reduction rate   : " << red_rate << "\n";
  cout << "output level control         : " << outlevel << "\n";
  cout << "\n";

  cout << "Results: \n\n";

  cout << "Steps (k):               " << k << "\n";
  cout << "Iterations (m):          " << m << "\n";
  cout << "Final state temperature: " << ck << "\n";
  cout << "\n";

  cout << "Last accepted solution: \n\n";

  print_solution(solution, Z);

  cout << "\n##################################################################"
          "###\n";
}
