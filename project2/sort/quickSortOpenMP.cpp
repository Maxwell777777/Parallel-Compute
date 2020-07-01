#include<iostream>
#include<omp.h>
#include<ctime>
#include <cstdlib>
#include<chrono>
#include<vector>
using namespace std;

int partition(vector<int> &arr, int low_index, int high_index)
{
    int left, right, middle;
    middle = arr[low_index];
    left= low_index + 1;
    right= high_index;
    while(left <= right)
    {
        while(left < high_index && middle >= arr[left])
            left++;
        while(middle < arr[right])
            right--;
        if(left < right)
        {
            int temp = arr[left];
            arr[left] = arr[right];
            arr[right] = temp;
        }
        else
        {
            int temp= arr[low_index];
            arr[low_index] = arr[right];
            arr[right]= temp;
            return(right);
        }
    }
    return(right);
}

void quicksort(vector<int> &arr, int low_index, int high_index)
{
    int right;

    if(low_index < high_index)
    {
        right = partition(arr, low_index, high_index);

        #pragma omp parallel sections
        {
            #pragma omp section
            {
                quicksort(arr, low_index, right - 1);
            }

            #pragma omp section
            {
                quicksort(arr, right + 1, high_index);
            }

        }
    }
}


int main(int argc,char** argv)
{
    int n = 1000000,i;
    int num_threads;
    n = atoi(argv[1]);
    num_threads = atoi(argv[2]);
    double start_time, run_time;
    srand(time(0));
    vector<int> arr;
    for(i=0;i < n;i++)
    {
        arr.push_back(rand() % n + 1);
    }
    omp_set_num_threads(num_threads);
    start_time = omp_get_wtime();
    quicksort(arr, 0, n - 1);
    run_time = omp_get_wtime() - start_time;
    cout<<"Elements of array after sorting \n";
    for(int i = 0; i < n; i++)
    {
        cout<<arr[i]<<"\t";
    }

    cout<<"\n";
    cout<<run_time;
    cout<<"\n";
}