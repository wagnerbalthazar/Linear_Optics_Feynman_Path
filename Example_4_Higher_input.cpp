/*************************************************************************
*
* Author:  Wagner F.Balthazar
* Date : Feb 13, 2023.
* Description : Simulation of boson sampling using Feynman's Integral path.

* Revision History :
1.0 - Feb 13, 2023 - First Edition
*
*
***************************************************************************/



/* Used libraries */

#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <complex>
#include <coroutine>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <windows.h>
#include <experimental/generator>
#include <valarray>
#include <chrono>

#include <bitset>

using namespace std;
double* all = new double[28];


double factorial(int n) {
    if (n == 0) {
        return 1;
    }
    int result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

int sum_vector(const vector<int>& v, int lower, int upper) {
    /* Function that sum all elements of a vector.

   Args:
   vector: The vector whose elements will be added.
   lower, upper: The integer whose given the limits of vector.

   Returns :
       int : The sum of the elements.
    */
    int sum = 0;
    for (int i = lower; i < upper; i++) {
        sum += v[i];
    }
    return sum;
}


int sum_vector_2(const vector<int>& v, int lower, int upper) {
    /* Function that sum all elements of a vector in step 2 by 2..

    Args:
    vector: The vector whose elements will be added.
    lower, upper: The integer whose given the limits of vector.

    Returns :
        int : The sum of the elements.
     */
    int sum = 0;
    for (int i = lower; i < upper; i += 2) {
        sum += v[i];
    }
    return sum;
}

complex <double> sum_vector_complex(const vector<complex<double>>& v, int lower, int upper) {

    /* Function that sum all elements of a vector with complex elements.

    Args:
    vector: The vector whose elements will be added.
    lower, upper: The integer whose given the limits of vector.

    Returns :
        int : The sum of the complex elements.
     */
    complex <double> sum = 0;
    for (int i = lower; i < upper; i++) {
        sum += v[i];
    }
    return sum;
}

vector<vector<complex<double>>> beamSplitter(double t1, double t2)
{
    /*Function to define the Beam Splitter(BS) matrix.

        Args:
        t1 : The argument related to the transmissivity of the beam splitter.
        t2: Relative phase between two arms of the Mach Zehnder interferometer.

        Returns :
        The 2 x 2 BS matrix.
     */

     //const double pi = acos(-1);
    const complex<double> i(0, 1);

    vector<vector<complex<double>>> BS = { {cos(t1),-(exp(-i * t2)) * sin(t1)},{sin(t1) * (exp(i * t2)),cos(t1)} };

    return BS;

}

bool is_vector_in_matrix(const vector<vector<int>>& matrix, const vector<int>& vec) {

    /*Function to verify if a vector is a row of a matrix.

        Args:
        t1 : The argument related to the transmissivity of the beam splitter.
        t2: Relative phase between two arms of the Mach Zehnder interferometer.

        Returns :
        True or false.
     */
    for (int i = 0; i < matrix.size(); i++) {
        bool match = true;
        for (int j = 0; j < vec.size(); j++) {
            if (matrix[i][j] != vec[j]) {
                match = false;
                break;
            }
        }
        if (match) {
            return true;
        }
    }
    return false;
}

int index_matrix(const vector<vector<int>>& matrix, const vector<int>& vec) {
    /*Return first index (row) of vector if it is a row in the matrix.
        Args:
        matrix: any matrix.
        vector: any vector.

        Returns :
        index or not found.
     */

    for (int i = 0; i < matrix.size(); i++) {
        if (matrix[i].size() == vec.size() && equal(begin(matrix[i]), end(matrix[i]), begin(vec))) {
            return i;
        }
    }
    return -1; // vector not found
}

bool count_negative(const vector<vector<int>>& matrix)
{
    /* Count the negative numbers. If one negative number is found the function stop or show that non-negative number is found.
        Args:
        matrix: any matrix.

        Returns :
        True (found) or false (there is no negative number).
     */

    for (const auto& i : matrix)
    {
        for (int j = 0; j < i.size(); j++)
        {
            if (i[j] < 0)
            {
                return true;
            }
        }
    }
    return false;
}

bool count_negative_aa(const vector<int>& vec)
{
    /* Count the negative numbers. If one negative number is found the function stop or show that non-negative number is found.
        Args:
        matrix: any vector.

        Returns :
        True (found) or false (there is no negative number).
     */

    for (int i : vec)
    {
        if (i < 0)
        {
            return true;
        }
    }

    return false;
}


vector<vector<int>> zerosMatrix(int n, int m) {

    /* Initialize matrix of zeros with n rows and m columns
        Args:
        n: rows.
        m: columns.

        Returns :
        matrix of zeros.
     */
    vector<vector<int>> result(n, vector<int>(m));

    for (vector<int>& row : result) {
        fill(row.begin(), row.end(), 0);
    }

    return result;
}


template <typename T>
int stack_vector(vector<vector<T>>& matrix, const vector<T>& vec, int dep) {
    /* Add a vector as a row in a matrix.
       Args:
       matrix: any matrix.
       vector: any vector.

       Returns:
       matrix with a new row.
    */
    int row = dep - 1;
    matrix.push_back(vec);
    return row;
}
void stack_vector_double(vector<vector<double>>& matrix, const vector<double>& vec) {
    /* Add a vector as a row in a matrix.
       Args:
       matrix: any matrix.
       vector: any vector.

       Returns:
       none.
    */
    matrix.push_back(vec);
}


vector<vector<complex<double>>> zerosMatrixComplex(int n, int m) {
    /* Initialize matrix of complex zeros with n rows and m columns
    Args:
    n: rows.
    m: columns.

    Returns :
    matrix of complex zeros.
 */
    vector<vector<complex<double>>> result(n, vector<complex<double>>(m));

    for (vector<complex<double>>& row : result) {
        fill(row.begin(), row.end(), 0);
    }

    return result;
}


// Function to create a vector of zeros
vector<int> zerosVector(int size)
{
    /* Initialize vector of zeros.
    Args:
    size:  any size.

    Returns :
    vector of zeros (int).
    */
    return vector<int>(size, 0);
}

vector<double> zerosVectorDouble(int size)
{
    /* Initialize matrix of zeros (double) with n rows and m columns
    Args:
    size:  any size.

    Returns :
    vector of zeros (double).
 */
    return vector<double>(size, 0);
}

vector<complex<double>> zerosVectorComplex(int size)
{
    /* Initialize matrix of zeros (complex).
    Args:
    size:  any size.

    Returns :
    vector of zeros (complex).
 */
    return vector<complex<double>>(size, complex<double>(0.0, 0.0));
}

vector<vector<int>> light_cone(vector<int> inp, vector<int> out, int dep)
{
    /// <summary>
    /// Determines what is the maximum occupation allowed in each interferometer arm. This is encoded into the matrix denoted as `general_matrix_loop`, 
    /// which is the used to construct the range list over which we are going to loop in our simulation procedure. This list is returned by this function, 
    /// and is denoted as `list_for_loop`.    
    /// </summary>
    /// <param name="inp"></param> The input (vector) of the interferometer.
    /// <param name="out"></param> The output (vector) of the interferometer.
    /// <param name="dep"></param> The depth of the interferometer.
    /// <returns></returns> The matrix of maximum (allowed) ranges in each green arm which we need to do interactions.

    int modes = inp.size();

    vector<int> first_layer_bs_occ(modes / 2, 0);
    int d = 0;
    for (int i = 0; i < modes; i += 2)
    {
        first_layer_bs_occ[d] = inp[i] + inp[i + 1];
        d += 1;
    }

    vector<vector<int>>  general_matrix_loop = zerosMatrix(modes, dep + 1);
    int t = 0;

    for (int i = 0; i < modes; i++) {
        general_matrix_loop[i][0] = inp[i];
        general_matrix_loop[i][dep] = out[i];
        if (2 * i < modes) {
            general_matrix_loop[i * 2][1] = first_layer_bs_occ[t];
            t += 1;
        }
    }

    for (int p = 2; p < dep - 1; p++) {
        if (p % 2 == 0) {
            for (int n = 1; n < (modes - 1); n += 2) {
                int start1 = max(n - p + 1, 0);
                int start2 = max(n - (dep - p), 0);
                int end1 = min(n + p + 1, modes);
                int end2 = min(n + dep - p, modes);
                int left_lcone = sum_vector(inp, start1, end1);
                int right_lcone = sum_vector(out, start2, end2);
                general_matrix_loop[n][p] = min(left_lcone, right_lcone);
            }

        }
        else {
            for (int n = 0; n < modes; n += 2) {
                int start1 = max(n - p + 1, 0);
                int start2 = max(n - (dep - p), 0);
                int end1 = min(n + p + 1, modes);
                int end2 = min(n + dep - p, modes);

                int left_lcone = sum_vector(inp, start1, end1);
                int right_lcone = sum_vector(out, start2, end2);
                general_matrix_loop[n][p] = min(left_lcone, right_lcone);
            }
        }

    }
    return general_matrix_loop;
}



vector<vector<int>> middle_light_cone_left(vector<int> inp, vector<int> out, int dep)
{
    /// <summary>
    /// /*Determines what is the maximum occupation allowed in each interferometer arm looking only for left of the circuit (input states). But in this function, we just consider part of the light cone from the horizontal position
    /// of the green wave guide (height) to bottom. This is encoded into the matrix denoted as `general_matrix_loop`, 
    /// which is the used to construct the range list over which we are going to loop in our simulation procedure. This list is returned by this function, 
    /// and is denoted as `list_for_loop`.
    /// <param name="inp"></param> The input (vector) of the interferometer.
    /// <param name="out"></param> The output (vector) of the interferometer.
    /// <param name="dep"></param> The depth of the interferometer.
    /// <returns></returns> The matrix of maximum (allowed) ranges in each green arm which we need to do interactions.
    /// 
    int modes = inp.size();

    vector<int> first_layer_bs_occ(modes / 2, 0);
    int d = 0;
    for (int i = 0; i < modes; i += 2)
    {
        first_layer_bs_occ[d] = inp[i] + inp[i + 1];
        d += 1;
    }

    vector<vector<int>>  general_matrix_loop = zerosMatrix(modes, dep + 1);
    int t = 0;
    for (int i = 0; i < modes; i++) {
        general_matrix_loop[i][0] = inp[i];
        general_matrix_loop[i][dep] = out[i];
        if (2 * i < modes) {
            general_matrix_loop[i * 2][1] = first_layer_bs_occ[t];
            t += 1;
        }
    }

    for (int p = 2; p < dep - 1; p++) {
        if (p % 2 == 0) {
            for (int n = 1; n < (modes - 1); n += 2) {
                int start = max(n - 1, 0);

                int end1 = min(n + p + 1, modes);
                int left_lcone = sum_vector(inp, start, end1);
                general_matrix_loop[n][p] = left_lcone;
            }

        }
        else {
            for (int n = 0; n < modes; n += 2) {
                int start = max(n, 0);
                int end1 = min(n + p + 1, modes);
                int left_lcone = sum_vector(inp, start, end1);
                general_matrix_loop[n][p] = left_lcone;
            }
        }
    }
    return general_matrix_loop;
}


vector<vector<int>> middle_light_cone_right(vector<int> inp, vector<int> out, int dep)
{
    /// <summary>
    /// /*Determines what is the maximum occupation allowed in each interferometer arm looking only for right of the circuit (output states). But in this function, we just consider part of the light cone from the horizontal position
    /// of the green wave guide (height) to bottom. This is encoded into the matrix denoted as `general_matrix_loop`, 
    /// which is the used to construct the range list over which we are going to loop in our simulation procedure. This list is returned by this function, 
    /// and is denoted as `list_for_loop`.
    /// <param name="inp"></param> The input (vector) of the interferometer.
    /// <param name="out"></param> The output (vector) of the interferometer.
    /// <param name="dep"></param> The depth of the interferometer.
    /// <returns></returns> The matrix of maximum (allowed) ranges in each green arm which we need to do interactions.
    /// 
    /// 
    int modes = inp.size();

    vector<int> first_layer_bs_occ(modes / 2, 0);
    int d = 0;
    for (int i = 0; i < modes; i += 2)
    {
        first_layer_bs_occ[d] = inp[i] + inp[i + 1];
        d += 1;
    }

    vector<vector<int>>  general_matrix_loop = zerosMatrix(modes, dep + 1);
    int t = 0;

    for (int i = 0; i < modes; i++) {
        general_matrix_loop[i][0] = inp[i];
        general_matrix_loop[i][dep] = out[i];
        if (2 * i < modes) {
            general_matrix_loop[i * 2][1] = first_layer_bs_occ[t];
            t += 1;
        }
    }

    for (int p = 2; p < dep - 1; p++) {
        if (p % 2 == 0) {
            for (int n = 1; n < (modes - 1); n += 2) {
                int start = max(n - 1, 0);
                int end2 = min(n + dep - p, modes);
                int right_lcone = sum_vector(out, start, end2);
                general_matrix_loop[n][p] = right_lcone;
            }

        }
        else {
            for (int n = 0; n < modes; n += 2) {
                int start = max(n, 0);
                int end2 = min(n + dep - p, modes);
                int right_lcone = sum_vector(out, start, end2);
                general_matrix_loop[n][p] = right_lcone;
            }
        }
    }
    return general_matrix_loop;
}



experimental::generator<vector<int>> odo_decr(vector<int> vector_odo) {

    /*  This function counts all possible occupations in one arm of BS given
        the total number of input photons in the same BS.If you consider a bunch of interferometers, max_ranges is
        a list which starts with all maximum photon occupation statesand decrease until that all outputs
        have no occupations.This functions works like an decreasing odometer.
        Args:
        vector: The maximum value of a loop (integer).

        Returns :
        vector: The bunch of list which corresponds of all possible number of occupattions.
    */


    vector<int> max_ranges = vector_odo;
    vector<int> item_for_loop(max_ranges.size(), 0);
    int pos = max_ranges.size() - 1;

    vector<int> greens;
    vector<int> greens_aux;

    greens_aux = max_ranges;

    int incr = 0;

    while (greens != item_for_loop) {
        if (incr == 0) {
            greens = max_ranges;
            incr += 1;
            co_yield greens;
        }
        else if (greens[pos] == 0) {
            greens[pos] = max_ranges[pos];
            pos -= 1;
            continue;
        }
        else {
            greens[pos] -= 1;
            pos = max_ranges.size() - 1;
            co_yield greens;
        }

    }
}

vector<double> generateRandomVector(int size) {

    /*  This function creates a vector with randow elements (double) between 0 and 1.
       Args:
       size: The number of elements of the vector.

       Returns :
       vector with ramdom elements.
   */

   // Create a random device
    random_device rd;

    // Seed a Mersenne Twister engine with the random device
    mt19937 engine(rd());

    // Create a uniform real distribution from 0 to 1
    uniform_real_distribution<double> dist(0, 1);

    // Create the vector
    vector<double> v(size);

    // Fill the vector with random numbers between 0 and 1
    for (int i = 0; i < size; i++) {
        v[i] = dist(engine);
    }

    return v;
}


vector<vector<vector<complex <double>>>> bs_parameters(bool random, int mod, int dep) {

    /* This function aims to allow us to choose the features of the BSs.Two possible options are allowed:
        50 / 50 (true) or random (false). For a personal setup, you must input the data by hand.

        //Args:
        //number_bs(int) : the number of BSs in the interferometer.
        //random(bool, optional) : Whether the beam splitters are 50 / 50 or random beam splitters. Defaults to False.

        Returns :
        Matrix (2 x 2 n), where n is the number of BSs, with the configuration of all BSS.
    */

    int number_bs = (mod / 2) * dep - int(dep / 2);

    vector<vector<vector<complex<double>>>> BSs;
    double pi = acos(-1);

    if (random) {
        vector<double> t11 = generateRandomVector(number_bs);
        vector<double> t22 = generateRandomVector(number_bs);
        vector<double> t1;
        vector<double> t2;
        for (double x : t11) {
            t1.push_back(x * pi / 2);
        }
        for (double y : t22) {
            t2.push_back(y * pi);
        }

        //std::cout << " t11 " << endl;
        //print_vector(t11);


        //std::cout << " t22 " << endl;
        //print_vector(t22);

        for (int i = 0; i < number_bs; i++) {
            // Get a matrix from the generateMatrix function
            //vector<vector<complex <double>>> matrix = beamSplitter(t1[i], t2[i]);
            BSs.push_back(beamSplitter(t1[i], t2[i]));
            /*BSs[0][2 * i] = matrix[0][0];
            BSs[0][2 * i + 1] = matrix[0][1];
            BSs[1][2 * i] = matrix[1][0];
            BSs[1][2 * i + 1] = matrix[1][1];*/
        }
    }
    else {

        double t1 = pi / 4;
        double t2 = 0;

        for (int i = 0; i < number_bs; i++) {
            // Get a matrix from the generateMatrix function
            //vector<vector<complex <double>>> matrix = beamSplitter(t1[i], t2[i]);
            BSs.push_back(beamSplitter(t1, t2));

        }


    }

    return BSs;

}


vector<vector<vector<complex <double>>>> bs_parameters_2(vector<double> par1, vector<double> par2, int mod, int dep) {

    /* This function aims to allow us to choose the features of the BSs.Two possible options are allowed:
        50 / 50 (true) or random (false). For a personal setup, you must input the data by hand.

        //Args:
        //number_bs(int) : the number of BSs in the interferometer.
        //random(bool, optional) : Whether the beam splitters are 50 / 50 or random beam splitters. Defaults to False.

        Returns :
        Matrix (2 x 2 n), where n is the number of BSs, with the configuration of all BSS.
    */

    int number_bs = (mod / 2) * dep - int(dep / 2);
    vector<vector<vector<complex<double>>>> BSs;
    double pi = acos(-1);
    vector<double> t1;
    vector<double> t2;
    for (double x : par1) {
        t1.push_back(x * pi / 2);
    }
    for (double y : par2) {
        t2.push_back(y * pi);
    }
    for (int i = 0; i < number_bs; i++) {
        // Get a matrix from the generateMatrix function
        //vector<vector<complex <double>>> matrix = beamSplitter(t1[i], t2[i]);
        BSs.push_back(beamSplitter(t1[i], t2[i]));
        /*BSs[0][2 * i] = matrix[0][0];
        BSs[0][2 * i + 1] = matrix[0][1];
        BSs[1][2 * i] = matrix[1][0];
        BSs[1][2 * i + 1] = matrix[1][1];*/
    }
    return BSs;
}



complex<double> amplitude_calculation_0(const vector<vector<int>>& new_matrix, int m, vector<vector<vector<complex<double>>>> bss, int dep)
{
    /*  This function computes the complex transition amplitudes for each possible occupation state
        considering only layer(less the first) of contraction. The function also finds the product of all
        amplitudes of this respective layer. It only works for the first layer of beam splitters.

        args:

        m: define what layer will be processed.
        new_matrix: the matrix that corresponding to the contraction layer with possible occupation states.
        bss: matrix with the arrays codifying the beam splitter matrices.
        ti: _description_
        tf: _description_
        amp: _description_
        pref: _description_
        namp1: _description_
        bsamp: _description_
        sumprodbsamp: the probability associated with a given inputand output combination.

        returns :
        probsamp: one complex amplitude transition for the layer.
     */


    int c = 0;
    auto pref = zerosVectorComplex((dep + 1) / 2);
    auto amp = zerosVectorComplex((dep + 1) / 2);
    auto bsamp = zerosVectorComplex((dep + 1) / 2);
    complex <double> namp1;
    int ti;
    int tf;
    int maxra;

    if (dep % 2 == 0) {
        maxra = dep - 1;

    }
    else {

        maxra = dep;
    }

    for (int ra = 0; ra < maxra; ra += 2) {
        if (m == 0) {


            int x1 = new_matrix[1][ra];
            int x2 = new_matrix[2][ra];
            int y1 = new_matrix[1][ra + 1];
            int y2 = new_matrix[2][ra + 1];

            if ((x1 + x2 != y1 + y2) or (y1 + y2 < x1) or (x1 + x2 < y1) or (x1 < 0) or (x2 < 0) or (y1 < 0) or (y2 < 0)) {

                amp[c] = 0;
                c += 1;
            }
            else {
                double n = 1;
                pref[c] = sqrt(factorial(x1) * factorial(x2) * factorial(y1) * factorial(y2));
                ti = max(max(0, y1 - x2), x1 - y2);
                tf = min(y1, x1);
                namp1 = (n / (factorial(ti) * factorial(y1 - ti) * factorial(x1 - ti) * factorial(x2 - y1 + ti))) * (pow(bss[c][0][0], ti)) * (pow(bss[c][0][1], (y1 - ti))) * (pow(bss[c][1][0], (x1 - ti))) * (pow(bss[c][1][1], (x2 - y1 + ti)));
                amp[c] = namp1;


                if (ti < tf) {
                    for (double i = ti; i < tf; i++) {
                        //double pm = (y1 - i);
                        namp1 = namp1 * (y1 - i) * (x1 - i) * bss[c][0][0] * bss[c][1][1] / ((i + 1) * (x2 - y1 + i + 1) * (bss[c][0][1]) * (bss[c][1][0]));
                        amp[c] += namp1;


                    }

                }
                c += 1;
            }


        }
    }

    for (int i = 0; i < amp.size(); i++) {
        bsamp[i] = pref[i] * amp[i];
    }

    complex <double> prodbsamp = 1;
    for (int j = 0; j < amp.size(); j++) {
        prodbsamp = prodbsamp * bsamp[j];
    }


    return prodbsamp;
}






complex<double> amplitude_calculation_n(const vector<vector<int>>& new_matrix, int m, vector<vector<vector<complex <double>>>> bss, int dep)
{
    /*  This function computes the complex transition amplitudes for each possible occupation state
    considering any layer (less the first and the last) of contraction. The function also finds the product of all
    amplitudes of these respective layer. The beam splitters are chosen one by one in the horizontal direction of its
    respective layer.
        args:

        m: define what layer will be processed.
        new_matrix : the matrix that corresponding to the contraction layer with possible occupation states.
        bss: matrix with the arrays codifying the beam splitter matrices.
        ti : _description_
        tf : _description_
        amp : _description_
        pref : _description_
        namp1 : _description_
        bsamp : _description_
        sumprodbsamp: the probability associated with a given inputand output combination.

        returns :
        probsamp: one complex amplitude transition for the layer.
     */

    vector<double> pref = zerosVectorDouble((dep));
    complex <double> namp1;
    vector<complex<double>> amp = zerosVectorComplex((dep));
    vector<complex<double>> bsamp = zerosVectorComplex((dep));

    int ti;
    int tf;

    int c = 0;
    int d = int(m / 2) * dep - int(dep / 2);


    for (int ra = 0; ra < dep; ra++) {
        if (ra % 2 == 0) {
            for (int rb = 1; rb < 2; rb++)
            {
                int x1 = new_matrix[rb][ra];
                int x2 = new_matrix[rb + 1][ra];
                int y1 = new_matrix[rb][ra + 1];
                int y2 = new_matrix[rb + 1][ra + 1];



                if ((x1 + x2 != y1 + y2) || (y1 + y2 < x1) || (x1 + x2 < y1) || (x1 < 0) || (x2 < 0) || (y1 < 0) || (y2 < 0)) {

                    amp[c] = 0;
                    c += 1;
                    d += 1;

                }
                else {

                    pref[c] = sqrt(factorial(x1) * factorial(x2) * factorial(y1) * factorial(y2));
                    ti = max(max(0, y1 - x2), x1 - y2);
                    tf = min(y1, x1);
                    double n = 1;
                    namp1 = n / (factorial(ti) * factorial(y1 - ti) * factorial(x1 - ti) * factorial(x2 - y1 + ti)) * pow(bss[d][0][0], ti) * pow(bss[d][0][1], (y1 - ti)) * pow(bss[d][1][0], (x1 - ti)) * pow(bss[d][1][1], (x2 - y1 + ti));

                    amp[c] = namp1;


                    if (ti < tf) {
                        for (int i = ti; i < tf; i++) {
                            namp1 = namp1 * double((y1 - i) * (x1 - i)) * bss[d][0][0] * bss[d][1][1] / (double((i + 1) * (x2 - y1 + i + 1)) * (bss[d][0][1]) * (bss[d][1][0]));
                            amp[c] += namp1;

                        }

                    }
                    c += 1;
                    d += 1;
                }
            }
        }
        else if (ra % 2 != 0)
        {
            for (int rb = 0; rb < 1; rb++)
            {
                int x1 = new_matrix[rb][ra];
                int x2 = new_matrix[rb + 1][ra];
                int y1 = new_matrix[rb][ra + 1];
                int y2 = new_matrix[rb + 1][ra + 1];

                if ((x1 + x2 != y1 + y2) || (y1 + y2 < x1) || (x1 + x2 < y1) || (x1 < 0) || (x2 < 0) || (y1 < 0) || (y2 < 0)) {

                    amp[c] = (0, 0);
                    c += 1;
                    d += 1;
                }
                else {
                    pref[c] = sqrt(factorial(x1) * factorial(x2) * factorial(y1) * factorial(y2));
                    ti = max(max(0, y1 - x2), x1 - y2);
                    tf = min(y1, x1);
                    double n = 1;
                    namp1 = n / (factorial(ti) * factorial(y1 - ti) * factorial(x1 - ti) * factorial(x2 - y1 + ti)) * pow(bss[d][0][0], ti) * pow(bss[d][0][1], (y1 - ti)) * pow(bss[d][1][0], (x1 - ti)) * pow(bss[d][1][1], (x2 - y1 + ti));

                    amp[c] = namp1;


                    if (ti < tf) {
                        for (int i = ti; i < tf; i++) {
                            //double pm = (y1 - i);
                            namp1 = namp1 * double((y1 - i) * (x1 - i)) * bss[d][0][0] * bss[d][1][1] / (double((i + 1) * (x2 - y1 + i + 1)) * (bss[d][0][1]) * (bss[d][1][0]));
                            amp[c] += namp1;

                        }

                    }
                    c += 1;
                    d += 1;
                }
            }
        }
    }

    for (int i = 0; i < amp.size(); i++) {
        bsamp[i] = pref[i] * amp[i];
    }

    complex <double> prodbsamp(1, 0);

    for (const auto& elem : bsamp) {
        prodbsamp *= elem;
    }

    return prodbsamp;
}

double probability_Feynman(vector<int> inp, vector<int> out, int dep, vector<vector<vector<complex <double>>>> bss)
{
    /// <summary>
    /// This function allows us to find out the probabilites of a boson sampling circuit for Clement's design using the Feynman approach.
    /// The code works for 'm' modes, if m >= 4 and for any depht >= 3.
    /// </summary>
    /// <param name="inp"></param> The vector with the input states.
    /// <param name="out"></param> The vector with the output states.
    /// <param name="dep"></param> The depth of the circuit.
    /// <param name="bss"></param> The beam splitters configurations.
    /// <returns></returns> A double with the probability.
    int sum_inp = sum_vector(inp, 0, inp.size());
    int sum_out = sum_vector(inp, 0, out.size());


    if (sum_inp != sum_out || inp.size() != out.size() || inp.size() % 2 != 0 || out.size() % 2 != 0) {
        throw std::runtime_error("The sum of the elements of vectors input and output must have the same number of photos, the same length and an even number os elements.");

    }
    else {

        int modes = inp.size();
        int number_bs = (modes / 2) * dep - static_cast<int>(dep / 2);
        vector<int> memory;

        if (dep % 2 != 0) {
            vector<vector<int>> general_matrix_loop_2 = light_cone(inp, out, dep);
            vector<vector<int>> general_matrix_loop_2_l = middle_light_cone_left(inp, out, dep);
            vector<vector<int>> general_matrix_loop_2_r = middle_light_cone_right(inp, out, dep);
            vector<vector<int>> pairs;
            vector<complex <double>> ampl;

            for (int m = 0; m < modes; m += 2) {
                vector<int> list_for_loop;
                vector<vector<int>>  matrix_contraction = zerosMatrix(3, dep + 1);


                if (m == 0) {

                    for (int n = 1; n < (dep - 1); n++) {
                        matrix_contraction[1][n] = general_matrix_loop_2[m][n];
                        matrix_contraction[2][n] = general_matrix_loop_2[m + 1][n];

                        if (n % 2 != 0) {
                            list_for_loop.push_back(matrix_contraction[1][n]);
                        }
                        else {
                            list_for_loop.push_back(matrix_contraction[2][n]);
                        }

                    }
                    for (auto item : odo_decr(list_for_loop)) {

                        vector<vector<int>>  matrix_contraction = zerosMatrix(3, dep + 1);
                        matrix_contraction[1][dep] = general_matrix_loop_2[m][dep];
                        matrix_contraction[2][dep] = general_matrix_loop_2[m + 1][dep];
                        matrix_contraction[1][0] = general_matrix_loop_2[m][0];
                        matrix_contraction[2][0] = general_matrix_loop_2[m + 1][0];
                        for (int n = 1; n < (dep); n++) {
                            if (n % 2 != 0) {
                                matrix_contraction[1][n] = item[n - 1];
                                matrix_contraction[1][n + 1] = matrix_contraction[1][n];
                                matrix_contraction[2][n] = matrix_contraction[1][n - 1] + matrix_contraction[2][n - 1] - matrix_contraction[1][n];

                            }
                            else if (n == dep - 1) {
                                matrix_contraction[2][n] = matrix_contraction[1][n + 1] + matrix_contraction[2][n + 1] - matrix_contraction[1][n];
                            }
                            else {
                                matrix_contraction[2][n] = item[n - 1];
                            }
                        }

                        if (count_negative(matrix_contraction) == 0) {

                            complex <double> amplit = amplitude_calculation_0(matrix_contraction, 0, bss, dep);
                            complex <double> n = 0;
                            if (abs(amplit) > pow(10, -15)) {
                                ampl.push_back(amplit);
                                vector <int> pairs2;

                                for (int i = 1; i < dep; i++) {
                                    pairs2.push_back(matrix_contraction[2][i]);
                                }
                                pairs.push_back(move(pairs2));

                            }
                        }
                    }
                }
                else if (m != 0 && m != modes - 2) {

                    vector<vector<int>> pairs3 = pairs;
                    vector<complex <double>> ampl2 = ampl;
                    pairs.clear();
                    ampl.clear();
                    matrix_contraction.clear();


                    vector<vector<int>>  matrix_contraction = zerosMatrix(3, dep + 1);
                    matrix_contraction[1][dep] = general_matrix_loop_2[m][dep];
                    matrix_contraction[2][dep] = general_matrix_loop_2[m + 1][dep];
                    matrix_contraction[1][0] = general_matrix_loop_2[m][0];
                    matrix_contraction[2][0] = general_matrix_loop_2[m + 1][0];

                    int c = 0;
                    int c_ii = -1;
                    int c_ii_aux = 0;

                    for (auto& ii : pairs3)
                    {

                        c_ii += 1;


                        vector<int> list_for_loop;
                        int len_ii = ii.size();

                        for (int t = 1; t < dep; t++)
                        {
                            matrix_contraction[0][t] = ii[t - 1];
                        }



                        for (int k = 1; k < (dep - 1); k++)
                        {
                            if (k == 1)
                            {
                                matrix_contraction[1][k] = general_matrix_loop_2_l[m][k];
                                list_for_loop.push_back(matrix_contraction[1][k]);

                            }
                            else if (k == 2)
                            {
                                matrix_contraction[2][k] = general_matrix_loop_2_l[m + 1][k];
                                list_for_loop.push_back(matrix_contraction[2][k]);

                            }
                            else if (k % 2 != 0 && k != 1) // && k != (dep - 2))
                            {
                                matrix_contraction[1][k] = min(general_matrix_loop_2_l[m][k] + sum_vector_2(ii, 0, k - 1),
                                    general_matrix_loop_2_r[m][k] + sum_vector_2(ii, k, len_ii));

                                list_for_loop.push_back(matrix_contraction[1][k]);

                            }
                            else if (k % 2 == 0 && k != 2)
                            {
                                matrix_contraction[2][k] = min(general_matrix_loop_2_l[m + 1][k] + sum_vector_2(ii, 0, k - 2),
                                    general_matrix_loop_2_r[m + 1][k] + sum_vector_2(ii, k + 1, len_ii));

                                list_for_loop.push_back(matrix_contraction[2][k]);

                            }
                            else if (k == dep - 2)
                            {
                                if (matrix_contraction[0][k + 1] + matrix_contraction[1][dep] + matrix_contraction[2][dep] - matrix_contraction[0][k] < 0)
                                {
                                    matrix_contraction[1][k] = min(general_matrix_loop_2_l[m][k] + sum_vector_2(ii, 0, k - 1),
                                        general_matrix_loop_2_r[m][k] + sum_vector_2(ii, k, len_ii));
                                    //matrix_contraction[1][k] = -1;

                                }
                                else
                                {
                                    matrix_contraction[1][k] = min(general_matrix_loop_2_l[m][k] + sum_vector_2(ii, 0, k - 1),
                                        general_matrix_loop_2_r[m][k] + sum_vector_2(ii, k, len_ii));
                                    /*matrix_contraction[1][k] = min(min(general_matrix_loop_2_l[m][k] + sum_vector_2(ii, 0, k - 1),
                                        general_matrix_loop_2_r[m][k] + sum_vector_2(ii, k, len_ii)),
                                        matrix_contraction[0][k + 1] + matrix_contraction[1][dep] + matrix_contraction[2][dep] - matrix_contraction[0][k]);*/

                                }

                                list_for_loop.push_back(matrix_contraction[1][k]);

                            }
                        }

                        vector<int> aa(int(pairs3[0].size() / 2) + 1, 0);

                        int w = 0;

                        for (int n = 1; n < dep; n += 2)
                        {
                            aa[w] = matrix_contraction[0][n] + (list_for_loop[n - 1]) - matrix_contraction[0][n + 1];

                            w += 1;
                        }

                        if (count_negative_aa(aa) == 0)
                        {

                            for (auto item : odo_decr(list_for_loop))
                            {
                                for (int n = 1; n < (dep - 1); n++) {
                                    if (n % 2 != 0) {
                                        matrix_contraction[1][n] = item[n - 1];
                                    }
                                    else {
                                        matrix_contraction[2][n] = item[n - 1];
                                    }

                                }
                                for (int t = 1; t < (dep); t++) {
                                    if (t % 2 != 0) {
                                        matrix_contraction[2][t] = matrix_contraction[1][t - 1] + matrix_contraction[2][t - 1] - matrix_contraction[1][t];

                                    }
                                    else {
                                        matrix_contraction[1][t] = matrix_contraction[1][t - 1] + matrix_contraction[0][t - 1] - matrix_contraction[0][t];
                                    }

                                }
                                matrix_contraction[2][dep - 1] = matrix_contraction[1][dep] + matrix_contraction[2][dep] - matrix_contraction[1][dep - 1];

                                if (count_negative(matrix_contraction) == 0)
                                {
                                    complex <double> ampl_prod = ampl2[c_ii] * amplitude_calculation_n(matrix_contraction, m, bss, dep);

                                    vector <int> pairs2;
                                    if (abs(ampl_prod) > pow(10, -15)) {
                                        for (int i = 1; i < dep; i++) {

                                            pairs2.push_back(matrix_contraction[2][i]);

                                        }

                                        if (is_vector_in_matrix(pairs, pairs2) == true)
                                        {
                                            int index = index_matrix(pairs, pairs2);
                                            ampl[index] += ampl_prod;

                                        }
                                        else
                                        {

                                            stack_vector(pairs, pairs2, dep);
                                            ampl.push_back(ampl_prod);

                                        }
                                    }

                                }

                            }
                        }
                    }
                }
                else if (m == modes - 2)
                {
                    vector<vector<int>> pairs3 = pairs;
                    vector<complex <double>> ampl2 = ampl;
                    pairs.clear();
                    ampl.clear();
                    vector<complex <double>> ampl3;

                    vector<vector<int>>  matrix_contraction = zerosMatrix(3, dep + 1);

                    matrix_contraction[1][dep] = general_matrix_loop_2[m][dep];
                    matrix_contraction[2][dep] = general_matrix_loop_2[m + 1][dep];
                    matrix_contraction[1][0] = general_matrix_loop_2[m][0];
                    matrix_contraction[2][0] = general_matrix_loop_2[m + 1][0];


                    int c = 0;
                    int c_ii = -1;
                    int c_ii_aux = 0;

                    for (int jj = 0; jj < pairs3.size(); jj++)
                    {
                        c_ii += 1;
                        vector<int> ii;
                        for (int mm = 0; mm < pairs3[0].size(); mm++) {

                            ii.push_back(pairs3[jj][mm]);

                        }

                        vector<int> list_for_loop;
                        int len_ii = ii.size();

                        for (int t = 1; t < dep; t++)
                        {
                            matrix_contraction[0][t] = ii[t - 1];
                        }

                        for (int k = 1; k < dep - 1; k++)
                        {
                            if (k == 1)
                            {
                                matrix_contraction[1][k] = general_matrix_loop_2_l[m][k];
                                list_for_loop.push_back(matrix_contraction[1][k]);


                            }
                            else if (k % 2 != 0 && k != 1)
                            {
                                matrix_contraction[1][k] = min(general_matrix_loop_2_l[m][k] + sum_vector_2(ii, 0, k - 1),
                                    general_matrix_loop_2_r[m][k] + sum_vector_2(ii, k, len_ii));

                                list_for_loop.push_back(matrix_contraction[1][k]);

                            }
                        }

                        vector<int> aa(static_cast<int>(pairs3[0].size() / 2), 0);

                        int w = 0;

                        for (int n = 1; n < dep; n += 2)
                        {
                            aa[w] = matrix_contraction[0][n] + (list_for_loop[w]) - matrix_contraction[0][n + 1];
                            w += 1;
                        }

                        if (count_negative_aa(aa) == 0)
                        {

                            vector<complex <double>> sum_ampl;
                            for (auto item : odo_decr(list_for_loop))
                            {

                                int c = 0;

                                for (int n = 1; n < (dep - 1); n += 2)
                                {
                                    matrix_contraction[1][n] = item[c];
                                    c += 1;

                                }
                                for (int t = 1; t < (dep); t++)
                                {
                                    if (t % 2 != 0) {
                                        matrix_contraction[2][t] = matrix_contraction[1][t - 1] + matrix_contraction[2][t - 1] - matrix_contraction[1][t];
                                        matrix_contraction[2][t + 1] = matrix_contraction[2][t];

                                    }
                                    else {
                                        matrix_contraction[1][t] = matrix_contraction[1][t - 1] + matrix_contraction[0][t - 1] - matrix_contraction[0][t];
                                    }

                                }

                                if (count_negative(matrix_contraction) == 0)
                                {
                                    complex <double> ampl_prod = ampl2[c_ii] * amplitude_calculation_n(matrix_contraction, m, bss, dep);
                                    if (abs(ampl_prod) > pow(10, -15)) {
                                        sum_ampl.push_back(ampl_prod);
                                    }

                                }

                            }

                            complex<double> a = sum_vector_complex(sum_ampl, 0, sum_ampl.size());
                            ampl3.push_back(a);

                        }

                    }
                    complex <double> amp_total = sum_vector_complex(ampl3, 0, ampl3.size());
                    complex<double> prob_comp = amp_total * conj(amp_total);
                    double prob = prob_comp.real();
                    return prob;
                }
            }
        }
        ///////////////////////////////// Even //////////////////////////////////////////////////////////
        else {
            vector<vector<int>> general_matrix_loop_2 = light_cone(inp, out, dep);
            vector<vector<int>> general_matrix_loop_2_l = middle_light_cone_left(inp, out, dep);
            vector<vector<int>> general_matrix_loop_2_r = middle_light_cone_right(inp, out, dep);

            vector<vector<int>> pairs;
            vector<complex <double>> ampl;

            for (int m = 0; m < modes; m += 2) {
                vector<int> list_for_loop;



                if (m == 0) {
                    vector<vector<int>>  matrix_contraction = zerosMatrix(3, dep);

                    for (int n = 1; n < (dep - 1); n++) {

                        if (n % 2 != 0) {
                            matrix_contraction[1][n] = general_matrix_loop_2[m][n];
                            list_for_loop.push_back(matrix_contraction[1][n]);
                        }
                        else {
                            matrix_contraction[2][n] = general_matrix_loop_2[m + 1][n];
                            list_for_loop.push_back(matrix_contraction[2][n]);
                        }

                    }

                    matrix_contraction[1][0] = general_matrix_loop_2[m][0];
                    matrix_contraction[2][0] = general_matrix_loop_2[m + 1][0];
                    matrix_contraction[1][dep - 1] = general_matrix_loop_2[m][dep];

                    for (auto item : odo_decr(list_for_loop)) {

                        for (int n = 1; n < (dep); n++) {
                            if (n % 2 != 0 && n != dep - 1) {
                                matrix_contraction[1][n] = item[n - 1];
                                matrix_contraction[1][n + 1] = matrix_contraction[1][n];
                                matrix_contraction[2][n] = matrix_contraction[1][n - 1] + matrix_contraction[2][n - 1] - matrix_contraction[1][n];

                            }
                            else if (n == dep - 1) {

                                matrix_contraction[2][n] = matrix_contraction[1][n - 1] + matrix_contraction[2][n - 1] - matrix_contraction[1][n];
                            }
                            else {
                                matrix_contraction[2][n] = item[n - 1];
                            }

                        }

                        if (count_negative(matrix_contraction) == 0) {
                            complex <double> amplit = amplitude_calculation_0(matrix_contraction, 0, bss, dep);

                            if (abs(amplit) > pow(10, -15)) {
                                ampl.push_back(amplit);
                                vector <int> pairs2;

                                for (int i = 1; i < dep; i++) {
                                    pairs2.push_back(matrix_contraction[2][i]);
                                }
                                pairs.push_back(move(pairs2));

                            }
                        }
                    }
                }
                else if (m != 0 && m != modes - 2) {

                    vector<vector<int>> pairs3 = pairs;
                    vector<complex <double>> ampl2 = ampl;
                    pairs.clear();
                    ampl.clear();
                    vector<vector<int>>  matrix_contraction = zerosMatrix(3, dep + 1);

                    matrix_contraction[1][0] = general_matrix_loop_2[m][0];
                    matrix_contraction[2][0] = general_matrix_loop_2[m + 1][0];
                    matrix_contraction[0][dep] = general_matrix_loop_2[m - 1][dep];
                    matrix_contraction[1][dep] = general_matrix_loop_2[m][dep];


                    int c = 0;
                    int c_ii = -1;
                    int c_ii_aux = 0;

                    for (auto& ii : pairs3)
                    {
                        c_ii += 1;

                        vector<int> list_for_loop;
                        int len_ii = ii.size();

                        for (int t = 1; t < dep; t++)
                        {
                            matrix_contraction[0][t] = ii[t - 1];
                        }

                        for (int k = 1; k < dep - 1; k++)
                        {
                            if (k == 1)
                            {
                                matrix_contraction[1][k] = general_matrix_loop_2_l[m][k];
                                list_for_loop.push_back(matrix_contraction[1][k]);


                            }
                            else if (k == 2)
                            {
                                matrix_contraction[2][k] = general_matrix_loop_2_l[m + 1][k];
                                list_for_loop.push_back(matrix_contraction[2][k]);

                            }
                            else if (k % 2 != 0 && k != 1)
                            {
                                matrix_contraction[1][k] = min(general_matrix_loop_2_l[m][k] + sum_vector_2(ii, 0, k - 1),
                                    general_matrix_loop_2_r[m][k] + sum_vector_2(ii, k, len_ii) + matrix_contraction[0][dep]);

                                list_for_loop.push_back(matrix_contraction[1][k]);

                            }
                            else if (k % 2 == 0 && k != 2)
                            {
                                matrix_contraction[2][k] = min(general_matrix_loop_2_l[m + 1][k] + sum_vector_2(ii, 0, k - 2),
                                    general_matrix_loop_2_r[m + 1][k] + sum_vector_2(ii, k + 1, len_ii) +
                                    matrix_contraction[0][dep]);

                                list_for_loop.push_back(matrix_contraction[2][k]);

                            }

                        }

                        vector<int> aa(static_cast<int>(pairs3[0].size() / 2 + 1), 0);

                        int w = 0;

                        for (int n = 1; n < dep; n += 2)
                        {
                            if (n < dep - 1) {
                                aa[w] = matrix_contraction[0][n] + (list_for_loop[n - 1]) - matrix_contraction[0][n + 1];
                                w += 1;
                            }
                            else {
                                aa[w] = matrix_contraction[0][dep] + matrix_contraction[1][dep] - matrix_contraction[0][dep - 1];
                            }
                        }

                        if (count_negative_aa(aa) == 0)
                        {

                            for (auto item : odo_decr(list_for_loop))
                            {
                                for (int n = 1; n < (dep - 1); n++) {
                                    if (n % 2 != 0) {
                                        matrix_contraction[1][n] = item[n - 1];
                                        matrix_contraction[2][n] = matrix_contraction[1][n - 1] + matrix_contraction[2][n - 1] - matrix_contraction[1][n];

                                    }
                                    else {
                                        matrix_contraction[2][n] = item[n - 1];
                                        matrix_contraction[1][n] = matrix_contraction[1][n - 1] + matrix_contraction[0][n - 1] - matrix_contraction[0][n];

                                    }

                                }
                                matrix_contraction[1][dep - 1] = matrix_contraction[0][dep] + matrix_contraction[1][dep] - matrix_contraction[0][dep - 1];
                                matrix_contraction[2][dep - 1] = matrix_contraction[2][dep - 2] + matrix_contraction[1][dep - 2] - matrix_contraction[1][dep - 1];


                                if (count_negative(matrix_contraction) == 0)
                                {
                                    complex <double> ampl_prod = ampl2[c_ii] * amplitude_calculation_n(matrix_contraction, m, bss, dep);
                                    if (abs(ampl_prod) > pow(10, -15)) {
                                        vector <int> pairs2;

                                        for (int i = 1; i < dep; i++) {

                                            pairs2.push_back(matrix_contraction[2][i]);

                                        }

                                        if (is_vector_in_matrix(pairs, pairs2) == true)
                                        {
                                            int index = index_matrix(pairs, pairs2);
                                            ampl[index] += ampl_prod;

                                        }
                                        else
                                        {

                                            stack_vector(pairs, pairs2, dep);
                                            ampl.push_back(ampl_prod);

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                else if (m == modes - 2) {
                    vector<vector<int>> pairs3 = pairs;
                    vector<complex <double>> ampl2 = ampl;
                    pairs.clear();
                    ampl.clear();
                    vector<complex <double>> ampl3;

                    vector<vector<int>>  matrix_contraction = zerosMatrix(3, dep + 1);

                    matrix_contraction[1][0] = general_matrix_loop_2[m][0];
                    matrix_contraction[2][0] = general_matrix_loop_2[m + 1][0];
                    matrix_contraction[0][dep] = general_matrix_loop_2[m - 1][dep];
                    matrix_contraction[1][dep] = general_matrix_loop_2[m][dep];

                    matrix_contraction[2][dep] = general_matrix_loop_2[m + 1][dep];
                    matrix_contraction[2][dep - 1] = matrix_contraction[2][dep];

                    int c = 0;
                    int c_ii = -1;
                    int c_ii_aux = 0;
                    for (int jj = 0; jj < pairs3.size(); jj++)
                    {
                        c_ii += 1;
                        vector<int> ii;
                        for (int mm = 0; mm < pairs3[0].size(); mm++) {

                            ii.push_back(pairs3[jj][mm]);

                        }

                        vector<int> list_for_loop;
                        int len_ii = ii.size();

                        for (int t = 1; t < dep; t++)
                        {
                            matrix_contraction[0][t] = ii[t - 1];
                        }
                        matrix_contraction[1][dep - 1] = matrix_contraction[0][dep] + matrix_contraction[1][dep] - matrix_contraction[0][dep - 1];

                        for (int k = 1; k < dep - 1; k += 2)
                        {
                            if (k == 1)
                            {
                                matrix_contraction[1][k] = general_matrix_loop_2_l[m][k];
                                list_for_loop.push_back(matrix_contraction[1][k]);


                            }
                            else if (k % 2 != 0 && k != 1)
                            {
                                matrix_contraction[1][k] = min(general_matrix_loop_2_l[m][k] + sum_vector_2(ii, 0, k - 1),
                                    general_matrix_loop_2_r[m][k] + sum_vector_2(ii, k, len_ii) + matrix_contraction[0][dep]);

                                list_for_loop.push_back(matrix_contraction[1][k]);

                            }

                        }

                        vector<int> aa(static_cast<int>(pairs3[0].size() / 2) + 1, 0);

                        int w = 0;

                        for (int n = 1; n < dep; n += 2)
                        {
                            if (n < dep - 2) {
                                aa[w] = matrix_contraction[0][n] + (list_for_loop[w]) - matrix_contraction[0][n + 1];
                                w += 1;

                            }
                            else {

                                aa[w] = matrix_contraction[1][dep - 1];
                            }

                        }

                        if (count_negative_aa(aa) == 0)
                        {
                            vector<complex <double>> sum_ampl;
                            for (auto item : odo_decr(list_for_loop))
                            {

                                int c = 0;

                                for (int n = 1; n < (dep - 2); n += 2)
                                {
                                    matrix_contraction[1][n] = item[c];
                                    c += 1;

                                }
                                for (int t = 1; t < (dep); t++)
                                {
                                    if (t % 2 != 0) {
                                        matrix_contraction[2][t] = matrix_contraction[1][t - 1] + matrix_contraction[2][t - 1] - matrix_contraction[1][t];
                                        matrix_contraction[2][t + 1] = matrix_contraction[2][t];

                                    }
                                    else {
                                        matrix_contraction[1][t] = matrix_contraction[1][t - 1] + matrix_contraction[0][t - 1] - matrix_contraction[0][t];
                                    }

                                }


                                if (count_negative(matrix_contraction) == 0)
                                {
                                    complex <double> ampl_prod = ampl2[c_ii] * amplitude_calculation_n(matrix_contraction, m, bss, dep);
                                    if (abs(ampl_prod) > pow(10, -15)) {
                                        sum_ampl.push_back(ampl_prod);
                                    }

                                }

                            }
                            complex<double> a = sum_vector_complex(sum_ampl, 0, sum_ampl.size());
                            ampl3.push_back(a);
                        }

                    }


                    complex <double> amp_total = sum_vector_complex(ampl3, 0, ampl3.size());
                    complex<double> prob_comp = amp_total * conj(amp_total);
                    double prob = prob_comp.real();
                    return prob;

                }

            }
        }
    }

}

int main(int argc, char* argv[]) {
    
    // Initial conditions
    int depth = 3;
    
    // Increasing the number of photons
    for (int j = 1; j < 21; ++j) {
        vector<int> input;
        for (int i = 0; i < 6; ++i) {
            input.push_back(j);
        }
        
        const vector<int> output = input;
        const int modes = input.size();

        vector<vector<vector<complex<double>>>> beam_splits = bs_parameters(true, modes, depth);

        // Probability and time Feyman.
        auto start = chrono::high_resolution_clock::now();
        cout << "Probabilty Feynman = " << probability_Feynman(input, output, depth, beam_splits) << endl; 
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
        double duration_seconds = static_cast<double>(duration.count()) / 1000000.0;
        cout << "Time taken Feynman: " << duration.count() / 1000000.0 << " seconds" << std::endl;

    }

    return 0;
}
