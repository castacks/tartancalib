#include <arpa/inet.h>

typedef size_t usize;
typedef int64_t i64;
typedef uint64_t u64;
typedef int32_t i32;
typedef uint32_t u32;
typedef int16_t i16;
typedef uint16_t u16;
typedef int8_t i8;
typedef uint8_t u8;

inline void write_one(const u8* data, FILE* file) {
fwrite(data, sizeof(u8), 1, file);
}

inline void write_one(const i8* data, FILE* file) {
fwrite(data, sizeof(i8), 1, file);
}

inline void write_one(const u16* data, FILE* file) {
u16 temp = htons(*data);
fwrite(&temp, sizeof(u16), 1, file);
}

inline void write_one(const i16* data, FILE* file) {
i16 temp = htons(*data);
fwrite(&temp, sizeof(i16), 1, file);
}

inline void write_one(const u32* data, FILE* file) {
u32 temp = htonl(*data);
fwrite(&temp, sizeof(u32), 1, file);
}

inline void write_one(const i32* data, FILE* file) {
i32 temp = htonl(*data);
fwrite(&temp, sizeof(i32), 1, file);
}

inline void write_one(const float* data, FILE* file) {
// TODO: Does this require a potential endian swap?
fwrite(data, sizeof(float), 1, file);
}
