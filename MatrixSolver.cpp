#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
 
using namespace std;
 
/** Prints the specified matrix onto the screen.
 */
void printMatrix(vector< vector<double> > matrix)
{
    for (int row = 0; row < matrix.size(); row++)
    {
        for (int column = 0; column < matrix[row].size(); column++)
            cout << left << setw(8) << matrix[row][column] << " ";
        cout << "\n";
    }
}
 
/** Multiplies each element in the given row by the scalar.
 */
void scale(vector<double> &row_val, vector<double> &row_unc, double scalar)
{
    for (int i = 0; i < row_val.size(); i++)
    {    
        row_val[i] *= scalar;
        row_unc[i] *= abs(scalar);
	}
}
 
/** Multiplies each element in the reference row by the scalar, and then
 * subtracts these values from the corresponding elements in the target row.
 */
void scalarSubtract (vector<double> &target_val, vector<double> &target_unc, double scalar, vector<double> ref_val, vector<double> ref_unc)
{
    for (int i = 0; i < target_val.size(); i++)
    {    
        target_val[i] -= scalar * ref_val[i];
        target_unc[i] = sqrt(target_unc[i] * target_unc[i] + scalar * scalar * ref_unc[i] * ref_unc[i]);
    }
}
 
/** Solves the given matrix.
 */
void solve(vector< vector<double> > &A, vector< vector<double> > &Unc)
{
    // Ensure that in A, all terms whose indices satisfy i = j are 
    // non-zero. This algorithm is not exhaustive and will likely  
    // fail if there are too many zeroes. However, its purpose is 
    // only to avert easily preventable div/0 errors.
    
    for (int i = 0; i < A.size(); i++)
    {
        // If i = j index element is zero, attempt to find replacement
        if (A[i][i] == 0) 
        {
            int column = 0;
            bool found = false;
            while (column < A.size() && !found)
            {
                if (A[i][column] != 0)
                {
                    if (A[column][i] != 0)
                    {
                        // If replacement row is valid, swap the rows in both
                        // matrix A and matrix B. Terminate search.
                        
                        swap(A[i], A[column]);
                        swap(Unc[i], Unc[column]);
                        found = true;
                    }
                }
                column++;
            }
        }
    }
    
    // GAUSS-JORDAN ELIMINATION:
     
    for (int i = 0; i < A.size(); i++)
    {
        // Multiply row i by scalar so that the i = j index element is 1.
         
        double scalar = 1/A[i][i];
        scale(A[i], Unc[i], scalar);
         
        // Scalar subtract row i from every row other than row i so that
        // row i contains the only non-zero element in column i.
         
        for (int row = 0; row < A.size(); row++)
        {
            if (row != i)
            {
                double scalar2 = A[row][i];
                scalarSubtract(A[row], Unc[row], scalar2, A[i], Unc[i]);
            }
        }
    }
}

/** Splits a string into a string array using the given delimiter.
 * Used to aid input parsing. Delimiter occurences are not included
 * in the returned array.
 */
vector<string> split(string str, string delim)
{
    // Declare Variables
    
    vector<string> params;
    string temp;
    size_t right = str.find(delim);
    
    // Iterate through string, finding occurrences of delimiter. 
    // Substring segment to array. Erase segment from source
    // string. Continue until no more delimiter occurrences.
    
    while (right != string::npos)
    {
        params.push_back(str.substr(0, right));
        str.erase(0, right + delim.size());
        right = str.find(delim);    
    }    
    params.push_back(str);    // Append the last remaining segment
    
    return params;
}


void inputMatrices(vector< vector<double> > &vals, vector< vector<double> > &uncs, string inpath, int rows, int cols)
{
	ifstream fin (inpath.c_str(), ios::in|ios::ate);
    if (fin.is_open())
    {       
        int size = fin.tellg();
        fin.seekg (0, ios::beg);
        vector<string>lines;    
        vector<string>params;    
        string line;
        while (getline(fin, line))        
            lines.push_back(line);               
        fin.close();
        
        int caret = 0;
        
        // Parse data
        
        vals = vector< vector<double> > (rows, vector<double>(cols));
    	uncs = vector< vector<double> > (rows, vector<double>(cols));
        
        while (caret < rows)
        {    
        	params = split(lines[caret], ",");
        	
        	for (int i = 0; i < params.size(); i++)
        		vals[caret][i] = atof(params[i].c_str());
        		
        	caret++;
        }
        
        caret++;
        
        while (caret < 2 * rows + 1)
        {
        	params = split(lines[caret], ",");
        	
        	for (int i = 0; i < params.size(); i++)
        		uncs[caret-rows-1][i] = atof(params[i].c_str());
        	
        	caret++;
        }       
	}
}

void outputMatrices(string outpath, vector< vector<double> > &vals, vector< vector<double> > &uncs)
{
	// Parse values into a stringstream
	
	 ostringstream data;
	 
	 for (int i = 0; i < vals.size(); i++)
	 {	 
	 	for (int j = 0; j < vals[i].size(); j++)
	 	{	 	
	 		data << vals[i][j];
	 		if (j < vals[i].size()-1)
	 			data << ",";
	 	}
	 	data << "\n";
	 }
	 
	 data << "\n";
	 
	 for (int i = 0; i < uncs.size(); i++)
	 {	 
	 	for (int j = 0; j < uncs[i].size(); j++)
	 	{	 	
	 		data << uncs[i][j];
	 		if (j < uncs[i].size()-1)
	 			data << ",";
	 	}
	 	data << "\n";
	 }
	 
	 // Write to file
	
	 ofstream fout(outpath.c_str(), ios_base::out | ios_base::trunc);

    // couldn't open it (disk error?); fail
    if (!fout.is_open())
    {    
        cout << "Could not write to " << outpath << endl;
    }

    fout.write(&(data.str()[0]), data.str().size());
    cout << "Wrote " << data.str().size() << " bytes to " << outpath << endl;
    fout.close();
}
 
int main()
{
    // Declaration of Variables
     
    vector< vector<double> > vals;
    vector< vector<double> > uncs;
    double temp;
    int rows, cols;
    string inpath, outpath;
     
    // Input    
    
    cout << "Rows = ";
    cin >> rows;
    cout << "Cols = ";
    cin >> cols;
    
    /*vals = vector< vector<double> > (rows, vector<double>(cols));
    uncs = vector< vector<double> > (cols, vector<double>(cols));
    
    cout << "Input data: " << endl;
    
    for (int i = 0; i < rows; i++)
    	for (int j = 0; j < cols; j++)
    		cin >> vals[i][j];
    		
	cout << "Input uncertainties: " << endl;
	
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			cin >> uncs[i][j];*/
			
	cout << "Input path = ";
	cin >> inpath;
	
	cout << "Output path = ";
	cin >> outpath;
	
	inputMatrices(vals, uncs, inpath, rows, cols);
			
	// Process
	
	solve(vals, uncs);
	
	// Output
	
	outputMatrices(outpath, vals, uncs);
         
    // End
    return 0;
}
