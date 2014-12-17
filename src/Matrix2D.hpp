#ifndef Matrix2D_hpp
#define Matrix2D_hpp


template<typename T>
class Matrix2D
{
    public:
      T*      data;
      size_t  nrows;
      size_t  ncols;   

    public:
      
      Matrix2D()
      {
          nrows = 0;
          ncols = 0;
          data  = 0;
      }
      
      Matrix2D(size_t m, size_t n)
      {
          allocate(m, n);
      }

      void allocate(size_t m, size_t n)
      {
          nrows = m;
          ncols = n;

          // make a column based allocation to be compatible with Fortran,
          // Matlab and especially BLAS routines that assume a column based
          // allocation of the matrix

          data = new T[ncols*nrows];

          if (data == NULL) 
          {
           printf("unable to allocate memory for Matrix (required: %ld bytes)\n", nrows*ncols*sizeof(double));
           return;
          }
      }

      void makeFrom(int m, int n, T* pdata)
      {
          nrows = m;
          ncols = n;

          memcpy(data, pdata, m*n);
      }

      ~Matrix2D()
      {
          delete[] data;
      } 

      // retrieve the j-th column
      T* operator[](size_t j)
      {
          return &data[nrows*j];
      }

      const T* operator[](size_t j) const
      {
          return &data[nrows*j];
      }

      T& operator()(size_t i, size_t j)
      {
          return data[nrows*j + i];
      }

     
      const T& operator()(size_t i, size_t j) const
      {
          return data[nrows*j + i];
      }

      void print() const
      {
          for (size_t i = 0; i < nrows; i++)
          {
              for (size_t j = 0; j < ncols; j++)
              {
                  cout << setw(25) << (*this)(i,j);
              }
              cout << endl;
          }
          cout << endl;
      }

      void saveToFile(const char* filename) const
      {
          fstream fout(filename, ios::out);
          fout.precision(16);
          fout.setf(ios::scientific, ios::floatfield);
          for (size_t i = 0; i < nrows; i++)
          {
              for (size_t j = 0; j < ncols; j++)
              {
                  fout << setw(25) << (*this)(i,j);
              }
              fout << endl;
          }
          fout.close();
      }

};

#endif

