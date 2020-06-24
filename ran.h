/****
	Random number generator from numerical recipies 3rd
****/

// basic type names (redefine if your bit lengths don't match)


typedef unsigned int Uint;

#ifdef _MSC_VER

typedef __int64 Llong; // 64 bit integer

typedef unsigned __int64 Ullong;

#else

typedef long long int Llong; // 64 bit integer

typedef unsigned long long int Ullong;

#endif

struct Ran {
Ullong u,v,w;
Ran(Ullong j) : v(4101842887655102017LL), w(1) {
u = j ^ v; int64();
v = u; int64();
w = v; int64();
}

  inline Ullong int64() {
    //returns an unsigned 64-bit integer
    //Ullmax = 2^(64)-1 = 18446744073709551615
    //1/Ullmax = 5.421010862427522e-20
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }
  inline double doub() {//returns a double-precision floating value in the range 0:0 to 1:0 
  return 5.42101086242752217E-20 * int64(); }

  inline double doub(double a,double b) {//returns a double-precision floating value in the range a to b
    return (b-a)*doub()+a;
  }
  inline Uint int32() { 
    //returns an unsigned 32-bit integer
    //Uimax=2^(32)-1
    return (Uint)int64(); }
};

