#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265

/*...................... Define Complex Number ...........................*/
typedef struct complex
{
    double real;
    double img;
}Complex;

/*...................... Evaluate Twiddle Base ...........................*/
void EvaluateTwiddleBase(int N,Complex TwiddleBase[]){

    int i;
    double val = PI /180;
    double theta;
    for(i = 0; i < N; i++)
    {
        theta = (2*180*i) / N;
        TwiddleBase[i].real = cos(theta*val);       // Cos takes radian
        TwiddleBase[i].img = -1 * sin(theta*val);   // Sin takes radian
    }

}

/*...................... Evaluate Twiddle Matrix Wn ...........................*/
void EvaluateTwiddleMatrix(int N,Complex TwiddleBase[],Complex TwiddleMatrix[N][N]){

    int i,j;
    // Using % as the values after N is Duplicated
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
        {
            TwiddleMatrix[i][j].real = TwiddleBase[(i*j)%N].real;
            TwiddleMatrix[i][j].img = TwiddleBase[(i*j)%N].img;
        }
    }
}

/*...................... Evaluate The Result of DFT ...........................*/
void EvaluateDFT(int N, Complex TwiddleMatrix[N][N],double X[], Complex DFT[]){

    int i,j;
    /* Assigning the result matrix with zeros */
        for(i = 0; i < N; i++)
        {
            DFT[i].real = 0;
            DFT[i].img = 0;
        }
    /* Evaluating the Result By Multiplying Wn[N][N] * X[N] */
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < N; j++)
            {
                DFT[i].real += TwiddleMatrix[i][j].real * X[j];
                DFT[i].img += TwiddleMatrix[i][j].img * X[j];
            }
        }
}

/*...................... Printing Discrete Fourier Transform X[K] ...........................*/
void PrintDFT(int N, Complex DFT[])
{   int i;

    /* Printing X[K] & Handling Real & Img if they equals Zeros */
    for(i = 0; i < N; i++)
    {   
        if(i == 0)
        printf("X[K} { ");

        if(i == N-1)
        {   
            if ( (int)(DFT[i].img * 100) == 0 )
                printf("%0.3lf, ",DFT[i].real);
            else if( (int)(DFT[i].real * 100) == 0 )
            {
                if(DFT[i].img > 0)
                printf("%0.3lfi ",DFT[i].img);
            else
                printf("%0.3lfi ",DFT[i].img);
            }
            else if(DFT[i].img > 0)
                printf("%0.3lf + %0.3lfi ",DFT[i].real,DFT[i].img);
            else
                printf("%0.3lf %0.3lfi ",DFT[i].real,DFT[i].img);
                printf("}");
        }
        else{
            if ( (int)(DFT[i].img * 100) == 0 )
                printf("%0.3lf, ",DFT[i].real);
            else if( (int)(DFT[i].real * 100) == 0 )
            {   
                if(DFT[i].img > 0)
                printf("%0.3lfi, ",DFT[i].img);
                else
                printf("%0.3lfi, ",DFT[i].img);
            }
            else if(DFT[i].img > 0)
                printf(" %0.3lf + %0.3lfi, ",DFT[i].real,DFT[i].img);
            else
                printf(" %0.3lf %0.3lfi, ",DFT[i].real,DFT[i].img);
        }
        
    }
    printf("\n");
}

/*...................... Printing Magnitude of X[K] ...........................*/
void PrintMagDFT(int N,double MagDFT[]){

        /* Printing Mag[X[K]] */
    int i;
    for(i = 0; i < N; i++)
    {   
        if(i == 0)
            printf("Mag[X[K]] { ");
        if(i == N-1)
        {   
            printf("%0.3lf ",MagDFT[i]);
            printf("}");
        }
        else{
            printf("%0.3lf ,",MagDFT[i]);
        }
    }
    printf("\n");
}

/*...................... Printing Phase of X[K] ...........................*/
void PrintPhaseDFT(int N,double PhaseDFT[])
{
    /* Printing Phase[X[K]] */
    int i;
    for(i = 0; i < N; i++)
    {   
        if(i == 0)
            printf("Phase[X[K]] { ");
        if(i == N-1)
        {   
            printf("%0.3lf ",PhaseDFT[i]);
            printf("}");
        }
        else{
            printf("%0.3lf, ",PhaseDFT[i]);
        }
        
    }
}

/*...................... Evaluate Magnitude of X[K] ...........................*/
void EvaluateMagintudeDFT(int N,Complex DFT[],double MagDFT[]){
    int i;

    for(i = 0; i < N; i++)
    {
        MagDFT[i] = sqrt( ( DFT[i].real * DFT[i].real ) + ( DFT[i].img * DFT[i].img ) );
    }
}

/*...................... Evaluate Phase of X[K] ...........................*/
void EvaluatePhaseDFT(int N,Complex DFT[],double PhaseDFT[]){
        int i;

    for(i = 0; i < N; i++)
    {      
        if ( (int)(DFT[i].img * 100) == 0)
        {
            if(DFT[i].real > 0)
                PhaseDFT[i] = 0;
            else
                PhaseDFT[i] = 180;
        }
        else if ( (int)(DFT[i].real * 100) == 0)
            {
                if(DFT[i].img > 0)
                    PhaseDFT[i] = 90;
                else
                    PhaseDFT[i] = -90;
            }
        else if(DFT[i].img > 0 && DFT[i].real > 0)
            PhaseDFT[i] = atan( (DFT[i].img) / (DFT[i].real) ) * 180 / PI;
        else if (DFT[i].img > 0 && DFT[i].real < 0)
            PhaseDFT[i] = 180 - ( atan( (DFT[i].img) / (-1 * (DFT[i].real) ) ) * 180 / PI );
        else if (DFT[i].img < 0 && DFT[i].real < 0)
            PhaseDFT[i] = -180 + ( atan ( (-1 * (DFT[i].img) ) / (-1 * (DFT[i].real) ) ) * 180 / PI  );
        else if (DFT[i].img < 0 && DFT[i].real > 0)
            PhaseDFT[i] = - ( atan( (-1 * (DFT[i].img) ) / (DFT[i].real) )  * 180 / PI );
    }
}

int main(int argc, char const *argv[])
{
    int N;    // Number of Points
    int i,j; // for loop iterators
    printf("How many N-Points do you want to calculate DFT? \n");
    scanf("%d", &N);    // Reading Number of Points
    printf("Enter the sequence form of X(N): Ex: 1 2 3 4 \n"); 
    double X[N],MagDFT[N],PhaseDFT[N];
    // Reading Input X[n] in Sequence Form
    for(i = 0; i < N; i++)
    {
        scanf("%lf", &X[i]);
    }
    Complex TwiddleBase[N];                                     // Twiddle Base array that includes n-1 elements
    Complex TwiddleMatrix[N][N];                                // Twiddle Matrix
    Complex DFT[N];
    EvaluateTwiddleBase(N, TwiddleBase);                        // Evaluate Twiddle Basic Elements
    EvaluateTwiddleMatrix(N, TwiddleBase, TwiddleMatrix);       // Evaluate Twiddle Matrix Wn
    EvaluateDFT(N, TwiddleMatrix, X, DFT);                      // Evaluate Discrete Fourier Transfor X[K]  
    EvaluateMagintudeDFT(N, DFT, MagDFT);                       // Evaluate Magnitude
    EvaluatePhaseDFT(N,DFT,PhaseDFT);                           // Evaluate Phase
    PrintDFT(N, DFT);                                           // Printing X[K]
    PrintMagDFT(N, MagDFT);                                     // Printing Magnitude 
    PrintPhaseDFT(N, PhaseDFT);                                 // Printing Phase

    return 0;
}
