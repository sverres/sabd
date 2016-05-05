# sabd


Case study on biodiversity data:

Given a forest with N treatment units. Each unit has a number of species of
interest. The owner decides to set aside Y units for preservation.

Which selection of Y units gives the largest number of species?

In this setup there is

- 155 units
- 133 species
- 10 units to be preserved


```cpp
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

float ck_0;             // initial temperature
int   Lk_0;             // initial step-length
float red_rate;         // temperature reduction rate
int   outlevel;         // output level control
```

Result of sample run showing 10 units with 101 species:

```
############## Simulated annealing results summary ##################

Parameters:

initial temperature,     ck_0: 35
iterations on each step, Lk_0: 50
temperature reduction rate   : 0.2
output level control         : 0

Results:

Steps (k):               54
Iterations (m):          2700
Final state temperature: 0.0002046102

Last accepted solution:

Z = 101    49  113  47  45  61  41  35  60  38  39      m: 2700

#####################################################################
```
