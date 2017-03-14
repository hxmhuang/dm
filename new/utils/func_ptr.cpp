
#include <vector>
#include <stdio.h>

typedef void (*func_ptr)(void** result, void* ptr,
			 void* alpha,  void* beta,
			 void* args,   int* n);

std::vector<func_ptr> ptrs;

extern "C"{
  
  void reg_func(void** ptr, int* id) {
    //printf("id=%d, %p\n", *id, *ptr);
    if(ptrs.size() < *id+1)
      ptrs.resize(*id+1);
    
    ptrs.at(*id) = (func_ptr)*ptr;
    //ptrs.insert(ptrs.begin() + *id, (func_ptr)*ptr);
  }

  void invoke(int* id, void** result,
	      void**ptr, void** alpha,
	      void** beta, void** args, int* n){
    
    if(*id==115){
      printf("result=%p\n", *result);
      printf("n=%d\n", *n);
    }
    
    *result = new int(1);
    //*(int*)(*result) = 1;
    
      //ptrs.at(*id)(result, *ptr, *alpha, *beta, *args, n);
  }
}
