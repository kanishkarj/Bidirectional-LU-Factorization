#include <bits/stdc++.h> 

using namespace std;

vector< double > b,y,x;
vector<double> sol;
vector<vector<double>> lower,upper;

void getRandomMatrix(int n,vector<vector<double>> &matrix) {
    srand(0);
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) {
            matrix[i][j] = rand()%150;
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
        for (int i=0 ; i<n ; i++) { 
        {
            y[i]=b[i];
            for (int j=0 ; j<i ; j++) { 
                y[i]-=lower[i][j]*y[j];
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
    
    // File 
    ofstream opfile;
    string filename = "doolittle-";
    ostringstream num_stream;
    num_stream << n;
    filename.append(num_stream.str());
    filename.append(".txt");
    opfile.open(filename);

    // The upper matrix
    upper.resize(n,vector<double>(n,0));
    // The Lower matrix
    lower.resize(n,vector<double>(n,0));
  

    double exec_time = clock();

    for(int k=0;k<n;k++)
    {
        for(int i=k+1;i<n;i++)
            matrix[i][k]=matrix[i][k]/matrix[k][k];

        for(int i=k+1;i<n;i++)
        {
            for(int j=k+1;j<n;j++)
            {
                matrix[i][j]=matrix[i][j]-matrix[i][k]*matrix[k][j];
            }
        }
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<=i;j++)
        {
            lower[i][j]=matrix[i][j];
            upper[j][i]=matrix[j][i];
        }
        lower[i][i]=1;
    }
    
    exec_time = clock() - exec_time;

    opfile<<endl;
    // Exectution time
    opfile<<"TIME TAKEN : "<<((double)exec_time/CLOCKS_PER_SEC)<<" milliseconds."<<endl;
    opfile<<endl;

    opfile<<"SIZE OF THE MATRIX "<<n<<" x "<<n<<endl<<endl;

    // Output
    opfile<<"Original Matrix : "<<endl;
    for (int i=0 ; i<n ; i++) { 
        for (int j=0 ; j<n ; j++) 
            opfile<<setw(10)<<fixed<<setprecision(0)<<matrix[i][j]<<"\t";
        opfile<<endl;  
    }

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
