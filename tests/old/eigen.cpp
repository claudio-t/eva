# include <eigen3/Eigen/Dense>
# include <eigen3/Eigen/Sparse>

# include <iostream>

/// Alias for type representing a real number
using real = double;

/// Alias for a type representing an index
using index_t = size_t;

/// Alias for a fixed size vector
template <int N, typename T = real>
using fixed_vector = Eigen::Matrix<T, N, 1>;

/// Alias for a fixed size matrix
template <int M, int N, typename T = real>
using fixed_matrix = Eigen::Matrix<T, M, N>;

/// Alias for a dense matrix used for computation
using dense_matrix = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>;

/// Alias for a dense vector used for computation
using dense_vector = Eigen::Matrix<real, Eigen::Dynamic, 1>;

/// Alias for a sparse matrix used for computation
using sparse_matrix = Eigen::SparseMatrix<real>; // Does not work with CRS!

/// Alias for a sparse vector used for computation (CCS format)
using sparse_vector = Eigen::SparseVector<real>;

int main(int argc, char * argv[])
{
    dense_vector  v1(5);
    v1 << 1, 2, 3, 4, 5;

    // sparse_vector v1;
    

    auto v2 = std::move(v1);

    std::cout << "v1 = " << v1 << std::endl;
    std::cout << "v2 = " << v2 << std::endl;

    std::cout << std::endl;
    
    std::cout << "v1 memaddr: " << v1.data() << std::endl;
    std::cout << "v2 memaddr: " << v2.data() << std::endl;

    
    return 0;
}
