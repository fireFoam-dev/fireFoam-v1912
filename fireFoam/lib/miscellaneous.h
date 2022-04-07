#ifndef MISCELLANEOUS_H
#define MISCELLANEOUS_H

#define DEBUG(x) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #x " = " << x << std::endl;
#define TRACE(s) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #s << std::endl; s;
#define DEBUGP(x) std::cout << "[p"<<Pstream::myProcNo()<<":"<< __FILE__ << ":" << __LINE__ << "] "<< #x " = " << x << std::endl;
#define TRACEP(s) std::cout << "[p"<<Pstream::myProcNo()<<":"<< __FILE__ << ":" << __LINE__ << "] "<< #s << std::endl; s;

#define CHECKIJ(index,jindex,variable) check(__FILE__,__LINE__,#variable,index,jindex,variable)
#define CHECKI(index,variable) check(__FILE__,__LINE__,#variable,index,variable)
#define CHECK(variable) check(__FILE__,__LINE__,#variable,variable)
#define CHECK0(variable) check0(__FILE__,__LINE__,#variable,variable)

#define COMM PETSC_COMM_WORLD

# ifdef LINUX
#    define ISNAN _isnan
#    define ISINF _isinf
# else
#    define ISNAN isnan
#    define ISINF isinf
# endif
# define ISZERO(x) fabs(x)<TINY
#define TINY 1e-16

#define VERBOSE

#endif
