#include <sys/mman.h>

//Aktually not Kilobytes: This would be Value*1000 -> but does not matter
#define Kilobytes(Value) ((Value)*1024LL)
#define Megabytes(Value) (Kilobytes(Value)*1024LL)
#define Gigabytes(Value) (Megabytes(Value)*1024LL)
#define Terabytes(Value) (Gigabytes(Value)*1024LL)

#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))

typedef size_t memory_index;


struct memory_arena
{
  memory_index Size;
  uint8* Base;
  memory_index Used;

};


#define PushStruct(Arena, type) (type *)PushSize_(Arena, sizeof(type))
#define PushArray(Arena, Count, type) (type *)PushSize_(Arena, (Count)*sizeof(type))

void * PushSize_(memory_arena *Arena, memory_index Size)
{
    Assert((Arena->Used + Size) <= Arena->Size);
    void *Result = Arena->Base + Arena->Used;
    Arena->Used += Size;
    
    return(Result);
}

memory_arena getMemoryArena(void *BaseAddress = (void *) 0, uint64 TotalSize = Megabytes(2))
{
    BaseAddress = (void *)(0);
    TotalSize = Megabytes(1);

    void *GameMemoryBlock = mmap(BaseAddress, (size_t)TotalSize,
                                 PROT_READ | PROT_WRITE,
                                 MAP_ANON | MAP_PRIVATE,
                                 -1, 0);

    Assert(GameMemoryBlock != MAP_FAILED);

    memory_arena Arena = {};
    Arena.Size = (size_t)TotalSize;
    Arena.Base = (uint8 *)GameMemoryBlock;

    return Arena;


}
