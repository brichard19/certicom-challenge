#ifndef _MANAGED_STACK_H
#define _MANAGED_STACK_H

#if !defined(__HIPCC__)
#include "hip_helper.h"
#endif

// Struct that can be passed to kernel
template<typename T> struct ManagedStack {
  T* ptr = nullptr;
  int* size = nullptr;
  int max_size = 0;
};

template<typename T> ManagedStack<T> stack_create(int max_size)
{
    ManagedStack<T> stack;

    stack.max_size = max_size;
    HIP_CALL(hipMallocManaged(&stack.ptr, sizeof(T) * max_size));
    HIP_CALL(hipMallocManaged(&stack.size, sizeof(int)));
    *stack.size = 0;

    return stack;
}

template<typename T> void stack_destroy(ManagedStack<T>& stack)
{
    HIP_IGNORE(hipFree(stack.size));
    HIP_IGNORE(hipFree(stack.ptr));
}

template<typename T> void stack_fill(ManagedStack<T>& stack, T* data, int count)
{
    memcpy(stack.ptr + *stack.size, data, sizeof(T) * count);
    *stack.size += count;
}

#if !defined(__HIPCC__)
template<typename T> __device__ T stack_pop(ManagedStack<T>& stack)
{
    int idx = *stack.size - 1;
    (*stack.size)--;

    return stack.ptr[idx];
}
#endif

// Device functions

#if defined(__HIPCC__)

template<typename T> __device__ __inline__ T stack_pop(ManagedStack<T>& stack)
{
    int idx = atomicSub(stack.size, 1) - 1;
    return stack.ptr[idx];
}

template<typename T> __device__ __inline__ void stack_push(ManagedStack<T>& stack, T& val)
{
    int idx = atomicAdd(stack.size, 1);
    stack.ptr[idx] = val;
}

#endif

#endif