/*************************************************************************
*
* Author:  Wagner F.Balthazar
* Date : Feb 13, 2023.
* Description : Additional code - Glynn's Formula.

* Revision History :
1.0 - Feb 22, 2023 - First Edition
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


using namespace std;



int factorial(int n) {
    /// <summary>
    ///   Function that computes the factorial of a given number.
    /// </summary>
    /// <param name="n"></param> The integer whose factorial is to be computed
    /// <returns></returns> The factorial of `number`.
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
    /// <summary>
    /// Function that sum all elements of a vector.
    /// </summary>
    /// <param name="v"></param>  The vector whose elements will be added
    /// <param name="lower"></param> The integer whose given lower limits of vector
    /// <param name="upper"></param> The integer whose given the upper of vector
    /// <returns></returns> The sum of the elements.

    int sum = 0;
    for (int i = lower; i < upper; i++) {
        sum += v[i];
    }
    return sum;
}

int sum_vector_2(const vector<int>& v, int lower, int upper) {
    /// <summary>
    /// Function that sum the elements of a vector in steps 2 by 2.
    /// </summary>
    /// <param name="v"></param>  The vector whose elements will be added
    /// <param name="lower"></param> The integer whose given lower limits of vector
    /// <param name="upper"></param> The integer whose given the upper of vector
    /// <returns></returns> The sum of the elements.
    int sum = 0;
    for (int i = lower; i < upper; i += 2) {
        sum += v[i];
    }
    return sum;
}

complex <double> sum_vector_complex(const vector<complex<double>>& v, int lower, int upper) {
    /// <summary>
    /// Function that sum all complex elements of a vector.
    /// </summary>
    /// <param name="v"></param>  The vector whose elements will be added
    /// <param name="lower"></param> The integer whose given lower limits of vector
    /// <param name="upper"></param> The integer whose given the upper of vector
    /// <returns></returns> The sum of the elements.
    /// 
    complex <double> sum = 0;
    for (int i = lower; i < upper; i++) {
        sum += v[i];
    }
    return sum;
}

void save_matrix_to_csv(vector<vector<double>> matrix, string filename) {
    /// <summary>
    /// Function to save a matrix in csv format
    /// </summary>
    /// <param name="matrix"></param> Matrix to save.
    /// <param name="filename"></param> Name of file.

    ofstream output_file(filename);
    for (auto row : matrix) {
        for (auto val : row) {
            output_file << val << ",";
        }
        output_file << endl;
    }
    output_file.close();
}

template <typename T>
void save_vector_to_csv(vector<T> vec, string filename) {
    /// <summary>
    /// Function to save a vector in csv format
    /// </summary>
    /// <param name="matrix"></param> Matrix to save.
    /// <param name="filename"></param> Name of file.
    ofstream output_file(filename);
    for (auto val : vec) {
        output_file << val << ",";
    }
    output_file << endl;
    output_file.close();
}

vector<vector<complex<double>>> beamSplitter(double t1, double t2) {
    /// <summary>
    /// *Function to define the Beam Splitter(BS) matrix.
    /// </summary>
    /// <param name="t1"></param> The argument related to the transmissivity of the beam splitter.
    /// <param name="t2"></param> Relative phase between two arms of the Mach Zehnder interferometer.
    /// <returns></returns> The 2 x 2 BS matrix.

    const complex<double> i(0, 1);

    vector<vector<complex<double>>> BS = { {cos(t1),-(exp(-i * t2)) * sin(t1)},{sin(t1) * (exp(i * t2)),cos(t1)} };

    return BS;

}

bool is_vector_in_matrix(const vector<vector<int>>& matrix, const vector<int>& vec) {
    /// <summary>
    /// Function to verify if a vector is a row of a matrix.
    /// </summary>
    /// <param name="matrix"></param> The matrix that you want to analyze.
    /// <param name="vec"></param> The vector that you want to know if it is in a matrix (row).
    /// <returns></returns> true or false.

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
    /// <summary>
    /// Return first index (row) of vector if it is a row in the matrix.
    /// </summary>
    /// <param name="matrix"></param> any matrix.
    /// <param name="vec"></param> any vector.
    /// <returns></returns> index or not found.

    for (int i = 0; i < matrix.size(); i++) {
        if (matrix[i].size() == vec.size() && equal(begin(matrix[i]), end(matrix[i]), begin(vec))) {
            return i;
        }
    }
    return -1; // vector not found
}

bool count_negative(const vector<vector<int>>& matrix) {
    /// <summary>
    /// Count the negative numbers. If one negative number is found the function stop or show that non-negative number is found.
    /// </summary>
    /// <param name="matrix"></param> any matrix.
    /// <returns></returns> True (found) or false (there is no negative number).
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

bool count_negative_aa(const vector<int>& vec) {
    /// <summary>
    /// Count the negative numbers. If one negative number is found the function stop or show that non-negative number is found.
    /// </summary>
    /// <param name="matrix"></param> any vector.
    /// <returns></returns> True (found) or false (there is no negative number).

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
    /// <summary>
    /// Initialize matrix of zeros with n rows and m columns (integers).
    /// </summary>
    /// <param name="n"></param> rows.
    /// <param name="m"></param> columns.
    /// <returns></returns> matrix of zeros.
    vector<vector<int>> result(n, vector<int>(m));

    for (vector<int>& row : result) {
        fill(row.begin(), row.end(), 0);
    }

    return result;
}

vector<vector<complex<double>>> zerosMatrixComplex(int n, int m) {
    /// <summary>
    /// Initialize matrix of zeros with n rows and m columns (double).
    /// </summary>
    /// <param name="n"></param> rows.
    /// <param name="m"></param> columns.
    /// <returns></returns> matrix of zeros.
    vector<vector<complex<double>>> result(n, vector<complex<double>>(m));

    for (vector<complex<double>>& row : result) {
        fill(row.begin(), row.end(), 0);
    }

    return result;
}

vector<int> zerosVector(int size) {
    /// <summary>
    /// Initialize vector of zeros with n rows and m columns (integers).
    /// </summary>
    /// <param name="n"></param> rows.
    /// <param name="m"></param> columns.
    /// <returns></returns> vector of zeros.
    return vector<int>(size, 0);
}

vector<double> zerosVectorDouble(int size)
{
    /// <summary>
    /// Initialize vector of zeros with n rows and m columns (double).
    /// </summary>
    /// <param name="n"></param> rows.
    /// <param name="m"></param> columns.
    /// <returns></returns> vector of zeros.
    return vector<double>(size, 0);
}

vector<complex<double>> zerosVectorComplex(int size) {
    /// <summary>
    /// Initialize vector of zeros with n rows and m columns (complex).
    /// </summary>
    /// <param name="n"></param> rows.
    /// <param name="m"></param> columns.
    /// <returns></returns> vector of zeros.
    return vector<complex<double>>(size, complex<double>(0.0, 0.0));
}

template <typename T>
int stack_vector(vector<vector<T>>& matrix, const vector<T>& vec, int dep) {
    /// <summary>
    /// Add a vector as a row in a matrix.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="matrix"></param> any matrix.
    /// <param name="vec"></param> any vector.
    /// <param name="dep"></param> depth of the circuit
    /// <returns></returns> matrix with a new row.
    int row = dep - 1;
    matrix.push_back(vec);
    return row;
}

vector<double> generate_random_vector(int size) {
    /// <summary>
    /// This function creates a vector with random uniform dsitribution elements (double) between 0 and 1.
    /// </summary>
    /// <param name="size"></param> The number of elements of the vector.
    /// <returns></returns> vector with ramdom elements.

    random_device rd;

    mt19937 engine(rd());

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

    /// <summary>
    /// This function aims to allow us to choose the features of the BSs.Two possible options are allowed:
    /// 50 / 50 (true) or random(false).For a personal setup, you must input the data by hand.
    /// </summary>
    /// <param name="random"></param> True - random beam splitters. False- 50/50 beam splitters
    /// <param name="mod"></param> The number of modes.
    /// <param name="dep"></param> The depth of the circuit.
    /// <returns></returns> Matrix(2 x 2n), where n is the number of BSs, with the configuration of all BSS.

    int number_bs = (mod / 2) * dep - int(dep / 2);

    vector<vector<vector<complex<double>>>> BSs;
    double pi = acos(-1);

    if (random) {
        vector<double> t11 = generate_random_vector(number_bs);
        vector<double> t22 = generate_random_vector(number_bs);
        vector<double> t1;
        vector<double> t2;
        for (double x : t11) {
            t1.push_back(x * pi / 2);
        }
        for (double y : t22) {
            t2.push_back(y * pi);
        }

        for (int i = 0; i < number_bs; i++) {
            BSs.push_back(beamSplitter(t1[i], t2[i]));
        }
    }
    else {

        double t1 = pi / 4;
        double t2 = 0;

        for (int i = 0; i < number_bs; i++) {
            BSs.push_back(beamSplitter(t1, t2));
        }
    }
    return BSs;
}

vector<vector<vector<complex <double>>>> bs_parameters_2(vector<double> par1, vector<double> par2, int mod, int dep) {

    /// <summary>
    /// This function aims to allow us to choose the features of the BSs by hand. For this, the beam splitters paramenters must be write in two vectors (par1, par2).
    /// </summary>
    /// <param name="par1"></param> A vector with the paramenters theta of the beam splitters.
    /// <param name="par2"></param> A vector with the paramenters phi of the beam splitters.
    /// <param name="mod"></param> The number of modes.
    /// <param name="dep"></param> The depth of the circuit.
    /// <returns></returns> Matrix(2 x 2n), where n is the number of BSs, with the configuration of all BSS.

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

int count_bits(unsigned int n) {
    /// <summary>
    /// This function returns the number of bits 1 in the binary vector.
    /// </summary>
    /// <param name="n"></param> the matrix that corresponding to the contraction layer with possible occupation states.
    /// <returns></returns> The number of 1s.
    ///     
    int count = 0;
    while (n) {
        count += n & 1;
        n >>= 1;
    }
    return count;
}

int find_one_position(int gray_code) {
    /// <summary>
    /// Find the position of 1 in a difference between two sequential gray codes
    /// </summary>
    /// <param name="gray_code"></param> The vector with gray code difference
    /// <returns></returns> The position of 1 in vector
    if (gray_code == 0) {
        return -1;
    }
    int position = 1;
    while ((gray_code & (position - 1)) != 0) {
        position <<= 1;
    }
    while (position > 0) {
        if ((gray_code & position) != 0) {
            return position / 2;
        }
        position /= 2;
    }
    return -1;
}

complex <double> perm_glynn(vector<vector<complex<double>>>& matrix) {
    /// <summary>
     /// This code calculates the permanent of a complex matrix using the Glynn's Formula.
     /// </summary>
     /// <param name="matrix"></param> any complex matrix.
     /// <returns></returns> The permanent of the matrix.
     /// 
    int n = matrix.size() - 1;
    complex<double> res(0, 0);

    for (int subs = 0; subs < (1 << n); ++subs) {
        vector<int> gray_vector;

        int gray = subs ^ (subs >> 1);
        std::bitset<32> b(gray);
        for (int i = 0; i < n; ++i) {

            gray_vector.push_back(b[n - i - 1]);

        }
        complex<double> prod(1, 0);
        for (int i = 0; i < n + 1; ++i) {
            complex<double> sum = matrix[0][i];
            for (int j = 1; j < n + 1; ++j) {
                if (gray_vector[j - 1]) {
                    sum -= matrix[j][i];
                }
                else {
                    sum += matrix[j][i];
                }

            }
            prod *= sum;
        }
        if ((n + 1) % 2 == 0) {
            if (sum_vector(gray_vector, 0, gray_vector.size()) % 2 == 0) {
                res -= prod;
            }
            else {
                res += prod;
            }
        }
        else {
            if (sum_vector(gray_vector, 0, gray_vector.size()) % 2 == 0) {
                res += prod;
            }
            else {
                res -= prod;
            }
        }
    }
    return res / pow(2, n);
}

vector<vector<complex<double>>> mat_unitary_circuit(vector<vector<vector<complex<double>>>> BSs, int dep, int mod) {

    /// <summary>
    /// This function builds the unitary matrix U which corresponds to the multi-mode interferometer.
    /// </summary>
    /// <param name="BSs"></param> Matrix with all beam splitter matrices of the circuit.
    /// <param name="dep"></param> Depth of the circuit.
    /// <param name="mod"></param> The number of modes of the circuit.
    /// <returns></returns>
    vector<vector<vector<complex<double>>>> layer1;
    for (int m = 0; m < dep; m++)
    {
        vector<vector<complex<double>>> layer = zerosMatrixComplex(mod, mod);
        if (m % 2 == 0) {
            int c = (m - m / 2);
            for (int n = 0; n < mod / 2; n++) {
                if (n == 0) {
                    layer[2 * n][2 * n] = BSs[c][0][0];
                    layer[2 * n][2 * n + 1] = BSs[c][0][1];
                    layer[2 * n + 1][2 * n] = BSs[c][1][0];
                    layer[2 * n + 1][2 * n + 1] = BSs[c][1][1];
                    c = m + (dep - dep / 2);
                }
                else {

                    layer[2 * n][2 * n] = BSs[c][0][0];
                    layer[2 * n][2 * n + 1] = BSs[c][0][1];
                    layer[2 * n + 1][2 * n] = BSs[c][1][0];
                    layer[2 * n + 1][2 * n + 1] = BSs[c][1][1];

                    c = c + dep;
                }
            }
        }
        else {
            int c = 0;
            c = m + (dep - dep / 2);
            for (int n = 1; n < mod / 2; n++) {

                layer[0][0] = 1;
                layer[2 * n - 1][2 * n - 1] = BSs[c][0][0];
                layer[2 * n - 1][2 * n] = BSs[c][0][1];
                layer[2 * n][2 * n - 1] = BSs[c][1][0];
                layer[2 * n][2 * n] = BSs[c][1][1];
                layer[mod - 1][mod - 1] = 1;

                c = c + dep;
            }
        }
        layer1.push_back(layer);
    }

    // Create a matrix to store the result of the multiplication
    vector<vector<complex<double>>> result = zerosMatrixComplex(mod, mod);
    vector<vector<complex<double>>> result2 = zerosMatrixComplex(mod, mod);
    vector<vector<complex<double>>> identity = zerosMatrixComplex(mod, mod);
    for (int i = 0; i < mod; i++) {
        result2[i][i] = complex<double>(1, 0);
    }

    // Perform matrix multiplication of all matrices in the vector
    for (int l = 0; l < dep; l++) {

        for (int i = 0; i < mod; i++) {
            for (int j = 0; j < mod; j++) {
                result[i][j] = (0, 0);
                for (int k = 0; k < mod; k++) {
                    result[i][j] += result2[i][k] * layer1[l][k][j];

                }
            }

        }
        result2 = result;
    }

    return result;

}

vector<vector<complex<double>>> Unitary_ST(vector<vector<complex<double>>> mat, vector<int> inp, vector<int> out) {
    /// <summary>
    /// This function defines a submatrix of U (unitary matrix f the circuit, see  mat_unitary_circuit), 
    /// which we call Ust, built considering the output (t) and take ti copies of the ith column, and input (S) taking 
    /// si copies of the ith row. 
    /// </summary>
    /// <param name="mat"></param> Unitary matrix U.
    /// <param name="inp"></param> Input state.
    /// <param name="out"></param> Output state.
    /// <returns></returns> The Ust matrix.
    /// 
    int mod = inp.size();
    vector<vector<complex<double>>> UT = zerosMatrixComplex(mod, sum_vector(out, 0, out.size()));
    int c = 0;
    for (int x = 0; x < out.size(); x++) {
        if (out[x] != 0) {
            vector<complex <double>> columnVector;
            for (int i = 0; i < mat.size(); i++) {
                columnVector.push_back(mat[i][x]);
            }
            int n = out[x];
            for (int j = 0; j < n; j++) {

                for (int i = 0; i < columnVector.size(); i++) {
                    UT[i][c] = columnVector[i];
                }
                c += 1;
            }
        }
    }
    vector<vector<complex<double>>> UST = zerosMatrixComplex(sum_vector(inp, 0, inp.size()), UT[0].size());
    int d = 0;
    for (int y = 0; y < inp.size(); y++) {
        if (inp[y] != 0) {
            vector<complex <double>> columnVector2;
            for (int i = 0; i < UT[0].size(); i++) {
                columnVector2.push_back(UT[y][i]);
            }
            int n1 = inp[y];
            for (int j = 0; j < n1; j++) {

                for (int i = 0; i < columnVector2.size(); i++) {
                    UST[d][i] = columnVector2[i];
                }
                d += 1;
            }
        }
    }
    return UST;

}

double probability_Glynn(vector<vector<complex<double>>> matrix, vector<int> inp, vector<int> out, int dep) {
    /// <summary>
    /// This function calculates the unitary matrix that corresponds to the multi-mode interfeformeter and the probality of transition given an input and an output
    /// using the Glynn's formula. 
    /// </summary>
    /// <param name="mat"><> /Matrix U which represents the interferometers given by clement's design.
    /// <param name="inp"><> Vector with input states (number of photons in each port).
    /// <param name="out"></param> Vector with output states (number of photons in each port).
    /// <param name="dep"></param> The depth of the circuit.
    /// <returns></returns> A double with the probability.
    int sum_inp = sum_vector(inp, 0, inp.size());
    int sum_out = sum_vector(out, 0, out.size());

    if (sum_inp != sum_out || inp.size() != out.size() || inp.size() % 2 != 0 || out.size() % 2 != 0) {
        throw std::runtime_error("The sum of the elements of vectors input and output must have the same number of photos, the same length and an even number os elements.");

    }
    else {
        int mod = inp.size();

        complex<double> permanent2 = perm_glynn(matrix);

        double prod = 1;
        for (int i = 0; i < mod; i++) {
            prod *= factorial(inp[i]) * factorial(out[i]);
        }
        complex<double> prob = permanent2 * conj(permanent2) / prod;

        return prob.real();
    }
}

int main(int argc, char* argv[]) {

    // Initial conditions
    int depth = 7;
    const vector<int> input = { 2,0,2,2,0,2 };
    const vector<int> output = { 1,3,3,1,0,0 };
    const int modes = input.size();
    vector<vector<vector<complex<double>>>> beam_splits = bs_parameters(true, modes, depth);
    vector<vector<complex<double>>> unitary_circuit = mat_unitary_circuit(beam_splits, depth, modes);
    vector<vector<complex<double>>> unitary_st = Unitary_ST(unitary_circuit, input, output);

    //Probability and time glynn.
    auto start2 = chrono::high_resolution_clock::now();
    cout << "probabilty glynn = " << probability_Glynn(unitary_st, input, output, depth) << endl;
    auto end2 = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration_cast<chrono::microseconds>(end2 - start2);
    double duration2_seconds = static_cast<double>(duration2.count()) / 1000000.0;
    cout << "time taken glynn: " << duration2.count() / 1000000.0 << " seconds" << std::endl;


    return 0;
}


