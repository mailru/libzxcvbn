#ifndef ZXCVBN_H
#define ZXCVBN_H

#include <stddef.h>
#include <stdint.h>
#include <sys/queue.h>
#include <time.h>

LIST_HEAD(zxcvbn_dict_head, zxcvbn_dict);

struct zxcvbn_node;

#define ZXCVBN_DATE_ONLY_YEAR   (1 << 0)
#define ZXCVBN_DATE_FULL_YEAR   (1 << 1)
#define ZXCVBN_DATE_SEPARATOR   (1 << 2)
/* date from passed list of dates */
#define ZXCVBN_DATE_FROM_LIST   (1 << 3)

#define zxcvbn_opts_init(o)     memset((o), 0, sizeof(*o))

typedef void *(*zxcvbn_malloc_t)(size_t size);
typedef void *(*zxcvbn_realloc_t)(void *ptr, size_t size);
typedef void (*zxcvbn_free_t)(void *ptr);

struct zxcvbn_opts {
    zxcvbn_malloc_t     malloc;
    zxcvbn_realloc_t    realloc;
    zxcvbn_free_t       free;
    const char         *symbols;
    unsigned int        max_matches_num;
};

struct zxcvbn_date {
    uint8_t     day;
    uint8_t     month;
    uint16_t    year;
    uint8_t     flags;
};

struct zxcvbn_spatial_graph {
    const char *data[256][9];
    double degree;
    unsigned int n_chars;
    unsigned int token_size;
    unsigned int n_coords;
};

struct zxcvbn_dict {
    struct zxcvbn *zxcvbn;
    LIST_ENTRY(zxcvbn_dict) list;
    int allocated;
    char name[PATH_MAX];
    struct zxcvbn_node *root;
};

struct zxcvbn {
    int allocated;
    void *(*zxcvbn_malloc)(size_t size);
    void *(*zxcvbn_realloc)(void *ptr, size_t size);
    void (*zxcvbn_free)(void *ptr);
    unsigned int n_symbols;
    char pack_table[256];
    unsigned int pack_table_size;
    struct zxcvbn_dict_head dict_head;
    struct zxcvbn_spatial_graph spatial_graph_qwerty;
    struct zxcvbn_spatial_graph spatial_graph_dvorak;
    struct zxcvbn_spatial_graph spatial_graph_keypad;
    struct zxcvbn_spatial_graph spatial_graph_macpad;
    unsigned int max_matches_num;
};

enum zxcvbn_match_type {
    ZXCVBN_MATCH_TYPE_DICT,
    ZXCVBN_MATCH_TYPE_SPATIAL,
    ZXCVBN_MATCH_TYPE_DIGITS,
    ZXCVBN_MATCH_TYPE_DATE,
    ZXCVBN_MATCH_TYPE_SEQUENCE,
    ZXCVBN_MATCH_TYPE_REPEAT,
    ZXCVBN_MATCH_TYPE_BRUTEFORCE,
};

#define ZXCVBN_MATCH_DESC_SEQ   (1 << 0)

struct zxcvbn_match {
    enum zxcvbn_match_type type;
    CIRCLEQ_ENTRY(zxcvbn_match) list;
    struct zxcvbn_spatial_graph *spatial_graph;
    union {
        struct zxcvbn_sequence     *seq;
        struct zxcvbn_date          date;
    };
    uint8_t                         flags;
    int i, j;
    unsigned int turns;
    unsigned int shifted;
    unsigned int rank;
    double entropy;
};

CIRCLEQ_HEAD(zxcvbn_match_head, zxcvbn_match);

struct zxcvbn_res {
    struct zxcvbn *zxcvbn;
    struct zxcvbn_match_head match_head;
    struct zxcvbn_match match_buf[32];
    struct zxcvbn_match *matches;
    unsigned int n_matches;
    unsigned int n_matches_reserved;
    double entropy;
};

struct zxcvbn *
zxcvbn_init(struct zxcvbn *zxcvbn_buf,
            void *(*zxcvbn_malloc)(size_t size),
            void *(*zxcvbn_realloc)(void *ptr, size_t size),
            void (*zxcvbn_free)(void *ptr),
            const char *symbols);

struct zxcvbn *
zxcvbn_init_ex(struct zxcvbn *zxcvbn, struct zxcvbn_opts *opts);

struct zxcvbn_dict *
zxcvbn_dict_init(struct zxcvbn *zxcvbn, struct zxcvbn_dict *dict_buf, const char *name);

int
zxcvbn_dict_add_word(struct zxcvbn_dict *dict, const char *word, unsigned int word_len, unsigned int rank);

void
zxcvbn_res_init(struct zxcvbn_res *res, struct zxcvbn *zxcvbn);

void
zxcvbn_res_release(struct zxcvbn_res *res);

int
zxcvbn_match(struct zxcvbn_res *res,
             const char *password, unsigned int password_len,
             char **words,         unsigned int words_num);

int
zxcvbn_match_ex(struct zxcvbn_res *res,
                const char *password,      unsigned int password_len,
                char **words,              unsigned int words_num,
                struct zxcvbn_date *dates, unsigned int dates_num);

const char *
zxcvbn_match_type_string(enum zxcvbn_match_type type);

void
zxcvbn_release(struct zxcvbn *zxcvbn);

#endif
