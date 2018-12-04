#include <bits/stdc++.h> 
#include <omp.h>
#include <ctime>
#include <ratio>
#include <chrono>

using namespace std;
using namespace std::chrono;

vector< double > b,y,x;
vector<double> sol;
vector<vector<double>> lower,upper;

void getRandomMatrix(int n,vector<vector<double>> &matrix) {
    srand(0);
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) {
            matrix[i][j] = rand()%199;
        }
    }
    sol.resize(n);
    for(int i=0 ; i<n ; i++){
        sol[i]=rand()%20;
    }
    b.resize(n);
    for (int i=0 ; i<n ; i++) { 

        for (int j=0 ; j<n ; j++) {
            b[i]+=matrix[i][j]*sol[j];
        }
    }
}
void forward_substitution(int n){

        y.resize(n);
        for (int i=0 ; i<n ; i++)
        {
            y[i]=b[i];
            for (int j=0 ; j<i ; j++) { 
                y[i]-=lower[i][j]*y[j];
              }
        }

}

void backward_substitution(int n){

    x.resize(n);

    for(int i=n-1;i>=0;i--)
        {
        x[i]=y[i];
            for(int j=i+1;j<n;j++)
                x[i]-=upper[i][j]*x[j];
            x[i]/=upper[i][i];
        }

}

void LU_Decomposition(vector<vector<double>>& matrix) 
{ 
    // Matrix Size
    int n = matrix.size();
    vector<vector<double>> temp;
    temp=matrix;
    // File 
    ofstream opfile;
    string filename = "parallel-";
    ostringstream num_stream;
    num_stream << n;
    filename.append(num_stream.str());
    filename.append(".txt");
    opfile.open(filename);

    // The upper matrix
    upper.resize(n,vector<double>(n,0));
    // The Lower matrix
    lower.resize(n,vector<double>(n,0));

    opfile<<endl;

    opfile<<"SIZE OF THE MATRIX "<<n<<" x "<<n<<endl<<endl;

    // Output
    opfile<<"Original Matrix : "<<endl;
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            opfile<<setw(10)<<fixed<<setprecision(0)<<matrix[i][j]<<"\t";
        opfile<<endl;  
    }
    
    auto start = high_resolution_clock::now();
 
    for(int k=0;k<n;k++)
    {
        #pragma omp parallel for
        for(int i=k+1;i<n;i++)
            matrix[i][k]=matrix[i][k]/matrix[k][k];

        #pragma omp parallel for
        for(int i=k+1;i<n;i++)
        {
            for(int j=k+1;j<n;j++)
            {
                matrix[i][j]=matrix[i][j]-matrix[i][k]*matrix[k][j];
            }
        }
    }


    #pragma omp parallel for
    for(int i=0;i<n;i++)
    {
        #pragma omp parallel for
        for(int j=0;j<=i;j++)
        {
            lower[i][j]=matrix[i][j];
            upper[j][i]=matrix[j][i];
        }
        lower[i][i]=1;
    }
    
    auto end = high_resolution_clock::now();
    auto time_span = duration_cast<duration<double>>(end - start);

    opfile<<endl;
    // Exectution time
    double time = time_span.count();
    opfile<<"TIME TAKEN : "<<time<<" seconds."<<endl;

    opfile<<endl;
    opfile<<"Lower Triangular : "<<endl; 
  
  
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            opfile<<setw(10)<<fixed<<setprecision(4)<<lower[i][j]<<"\t";
        opfile<<endl;  
    }

    opfile<<endl<<endl;

    opfile<<"Upper Triangular : "<<endl; 
  
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            opfile<<setw(10)<<fixed<<setprecision(4)<<upper[i][j]<<"\t"; 
        opfile<<endl;  
    } 

    // Check if the LU decomposition is correct by multiplying it back.
    for (int i=0 ; i<n ; i++) 
    { 
        for (int j = 0 ; j < n ; j++) 
        { 
            matrix[i][j] = 0; 
            for (int k = 0 ; k < n ; k++) 
                matrix[i][j] += lower[i][k] * upper[k][j]; 
        } 
    } 

    opfile<<endl<<endl;
    opfile<<"Original matrix after multiplication : "<<endl; 
    
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            opfile<<setw(10)<<fixed<<setprecision(0)<<matrix[i][j]<<"\t"; 
        opfile<<endl;  
    } 

    forward_substitution(n);
    backward_substitution(n);

    opfile<<"For B = : ";
    opfile<<endl;
    for (int j=0 ; j<b.size() ; j++) 
        opfile<<setw(10)<<fixed<<setprecision(0)<<b[j]<<"\t"; 
    opfile<<endl;

    opfile<<"X = : ";
    opfile<<endl;
    for (int j=0 ; j<x.size() ; j++) 
        opfile<<setw(10)<<fixed<<setprecision(0)<<x[j]<<"\t"; 
    opfile<<endl;

    cout<<"Data written to "<<filename;
}

int main(int argc, char** args) 
{ 
    int n = atoi(args[1]);
    // cout<<n;
    vector<vector<double>> matrix(n,vector<double>(n,0));
    getRandomMatrix(n,matrix);

    LU_Decomposition(matrix); 
    
    return 0; 
} 
