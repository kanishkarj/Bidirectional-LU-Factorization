#include <bits/stdc++.h> 

using namespace std; 

void LU_Decomposition(vector<vector<int>> matrix) 
{ 
    // Matrix Size
    int n = matrix.size();

    // The upper matrix
    vector<vector<int>> lower(n,vector<int>(n,0));
    // The Lower matrix
    vector<vector<int>> upper(n,vector<int>(n,0));
  
    // Decomposing matrix into Upper and Lower 
    // triangular matrix 
    for (int i=0 ; i<n; i++) { 
  
        // Upper Triangular 
        for (int k=i; k<n ; k++) { 
  
            // Summation of L(i, j) * U(j, k) 
            int sum = 0; 
            for (int j=0 ; j<i ; j++) 
                sum += (lower[i][j] * upper[j][k]); 
  
            // Evaluating U(i, k) 
            upper[i][k] = matrix[i][k] - sum; 
        } 
  
        // Lower Triangular 
        for (int k=i ; k<n ; k++) { 
            if (i == k) 
                lower[i][i] = 1; // Diagonal as 1 
            else { 
  
                // Summation of L(k, j) * U(j, i) 
                int sum = 0; 
                for (int j=0 ; j<i ; j++) 
                    sum += (lower[k][j] * upper[j][i]); 
  
                // Evaluating L(k, i) 
                lower[k][i] = (matrix[k][i] - sum) / upper[i][i]; 
            } 
        } 
    } 
  
    // Output
    cout<<"Lower Triangular : "<<endl; 
  
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            cout<<lower[i][j]<<"\t";
        cout<<endl;  
    }

    cout<<endl<<endl;

    cout<<"Upper Triangular : "<<endl; 
  
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            cout<<upper[i][j]<<"\t"; 
        cout<<endl;  
    } 
} 

// Driver code 
int main() 
{ 
    vector<vector<int>> matrix = 
    { 
        { 2, -1, -2 }, 
        { -4, 6, 3 }, 
        { -4, -2, 8 } 
    }; 
  
    LU_Decomposition(matrix); 
    return 0; 
} 