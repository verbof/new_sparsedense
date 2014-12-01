#ifndef Matrix2D_hpp
#define Matrix2D_hpp


template<typename T>
class Matrix2D
{
    public:
      T**     data;
      size_t  rows;
      size_t  cols;   

    public:
      
      Matrix2D()
      {
          rows = 0;
          cols = 0;
          data = 0;
      }
      
      Matrix2D(size_t n, size_t m)
      {
          allocate(n, m);
      }

      void allocate(size_t n, size_t m)
      {
          rows = n;
          cols = m;

          // make a column based allocation to be compatible with Fortran,
          // Matlab and especially BLAS routines that assume a column based
          // allocation of the matrix

          data    = new T*[cols];
          data[0] = new T[cols*rows];

          if (data == NULL) 
          {
           printf("unable to allocate memory for Matrix (required: %ld bytes)\n", rows*cols*sizeof(double));
           return EXIT_FAILURE;
          }


          for (size_t i = 1; i < cols; i++)
          {
              data[i] = data[i-1] + rows;
          }
      }

      ~Matrix2D()
      {
          delete[] data[0];
          delete[] data;
      } 

      // retrieve the j-th column
      T* operator[](size_t j)
      {
          return data[j];
      }

      const T* operator[](size_t j) const
      {
          return data[j];
      }

      T& operator()(size_t i, size_t j)
      {
          return data[j][i];
      }

     
      const T& operator()(size_t i, size_t j) const
      {
          return data[j][i];
      }

      void print() const
      {
          for (size_t i = 0; i < rows; i++)
          {
              for (size_t j = 0; j < cols; j++)
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
          for (size_t i = 0; i < rows; i++)
          {
              for (size_t j = 0; j < cols; j++)
              {
                  fout << setw(25) << (*this)(i,j);
              }
              fout << endl;
          }
          fout.close();
      }

};

#endif

