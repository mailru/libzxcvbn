#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "zxcvbn.h"

#ifndef ZXCVBN_PASSWORD_LEN_MAX
#define ZXCVBN_PASSWORD_LEN_MAX 256
#endif

struct zxcvbn_dict;

// prefix tree
struct zxcvbn_node {
    struct zxcvbn_node **children;
    int rank;
};

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

static inline void *
__malloc(struct zxcvbn *zxcvbn, size_t size)
{
    return (*zxcvbn->zxcvbn_malloc)(size);
}

static inline void *
__realloc(struct zxcvbn *zxcvbn, void *ptr, size_t size)
{
    return (*zxcvbn->zxcvbn_realloc)(ptr, size);
}

static inline void
__free(struct zxcvbn *zxcvbn, void *ptr)
{
    return (*zxcvbn->zxcvbn_free)(ptr);
}

const char *
zxcvbn_match_type_string(enum zxcvbn_match_type type)
{
    switch (type) {
    case ZXCVBN_MATCH_TYPE_DICT:
        return "dict";
    case ZXCVBN_MATCH_TYPE_SPATIAL:
        return "spatial";
    case ZXCVBN_MATCH_TYPE_DIGITS:
        return "digits";
    case ZXCVBN_MATCH_TYPE_DATE:
        return "date";
    case ZXCVBN_MATCH_TYPE_BRUTEFORCE:
        return "bruteforce";
    default:
        assert(0);
    }
}

static unsigned int
get_align_coords(int coords[2][8], int x, int y)
{
    coords[0][0] = x - 1;
    coords[1][0] = y;

    coords[0][1] = x - 1;
    coords[1][1] = y - 1;

    coords[0][2] = x;
    coords[1][2] = y - 1;

    coords[0][3] = x + 1;
    coords[1][3] = y - 1;

    coords[0][4] = x + 1;
    coords[1][4] = y;

    coords[0][5] = x + 1;
    coords[1][5] = y + 1;

    coords[0][6] = x;
    coords[1][6] = y + 1;

    coords[0][7] = x - 1;
    coords[1][7] = y + 1;

    return 8;
}

static unsigned int
get_liner_coords(int coords[2][8], int x, int y)
{
    coords[0][0] = x - 1;
    coords[1][0] = y;

    coords[0][1] = x + 1;
    coords[1][1] = y;

    return 2;
}

static unsigned int
get_slant_coords(int coords[2][8], int x, int y)
{
    coords[0][0] = x - 1;
    coords[1][0] = y;

    coords[0][1] = x;
    coords[1][1] = y - 1;

    coords[0][2] = x + 1;
    coords[1][2] = y - 1;

    coords[0][3] = x + 1;
    coords[1][3] = y;

    coords[0][4] = x;
    coords[1][4] = y + 1;

    coords[0][5] = x - 1;
    coords[1][5] = y + 1;

    return 6;
}

static void
make_spatial_graph_iter(struct zxcvbn_spatial_graph *spatial_graph,
                        const char **kb, unsigned int x_size, unsigned int y_size,
                        unsigned int token_size, unsigned int n_coords,
                        unsigned int get_coords(int coords[2][8], int x, int y))
{
#define  KB_XY(x, y) \
    kb[x_size *y + x]

    int x, y, i, coords[2][8], n_coords_get;
    const char *s;

    assert(y_size >= 2);

    spatial_graph->token_size = token_size;
    spatial_graph->n_coords = n_coords + 1;

    for (y = 1; y < y_size - 1; ++y) {
        for (x = 1; x < x_size - 1; ++x) {
            if ((s = KB_XY(x, y)) == NULL)
                continue;

            n_coords_get = get_coords(coords, x, y);
            assert(n_coords_get == n_coords);
            for (; *s != '\0'; ++s) {
                ++spatial_graph->n_chars;
                for (i = 0; i < n_coords; ++i)
                    if ((spatial_graph->data[*s][i] = KB_XY(coords[0][i], coords[1][i])) != NULL)
                        ++spatial_graph->degree;
                spatial_graph->data[*s][n_coords] = KB_XY(x, y);
            } 
        }
    }

    for (y = 0; y < 256; ++y) {
        for (x = 0; x < spatial_graph->n_coords; ++x) {
            if (spatial_graph->data[y][x] == NULL)
                spatial_graph->data[y][x] = "\xff\xff";
        }
    }

    spatial_graph->degree /= spatial_graph->n_chars;

#undef KB_XY
}

static int
make_spatial_graph(struct zxcvbn *zxcvbn)
{
#define KEYBRD_X_SIZE 16
#define KEYBRD_Y_SIZE 6
#define KEYPAD_X_SIZE 6
#define KEYPAD_Y_SIZE 7
#define ALPHABET_X_SIZE 28
#define ALPHABET_Y_SIZE 3

    static const char *qwerty[KEYBRD_Y_SIZE][KEYBRD_X_SIZE] = {
        { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL,  NULL },
        { NULL, "`~", "1!", "2@", "3#", "4$", "5%", "6^", "7&", "8*", "9(", "0)", "-_",  "=+", NULL,  NULL },
        { NULL, NULL, "qQ", "wW", "eE", "rR", "tT", "yY", "uU", "iI", "oO", "pP", "[{",  "]}", "\\|", NULL },
        { NULL, NULL, "aA", "sS", "dD", "fF", "gG", "hH", "jJ", "kK", "lL", ";:", "'\"", NULL, NULL,  NULL }, 
        { NULL, NULL, "zZ", "xX", "cC", "vV", "bB", "nN", "mM", ",<", ".>", "/?", NULL,  NULL, NULL,  NULL },
        { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL,  NULL }
    };

    static const char *dvorak[KEYBRD_Y_SIZE][KEYBRD_X_SIZE] = {
        { NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, },
        { NULL, "`~", "1!",  "2@", "3#", "4$", "5%", "6^", "7&", "8*", "9(", "0)", "[{", "]}", NULL,  NULL, },
        { NULL, NULL, "'\"", ",<", ".>", "pP", "yY", "fF", "gG", "cC", "rR", "lL", "/?", "=+", "\\|", NULL, },
        { NULL, NULL, "aA",  "oO", "eE", "uU", "iI", "dD", "hH", "tT", "nN", "sS", "-_", NULL, NULL,  NULL, },
        { NULL, NULL, ";:",  "qQ", "jJ", "kK", "xX", "bB", "mM", "wW", "vV", "zZ", NULL, NULL, NULL,  NULL, },
        { NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, },
    };

    static const char *keypad[KEYPAD_Y_SIZE][KEYPAD_X_SIZE] = {
        { NULL, NULL, NULL, NULL, NULL, NULL, },
        { NULL, NULL, "/",  "*",  "-",  NULL, },
        { NULL, "7",  "8",  "9",  "+",  NULL, },
        { NULL, "4",  "5",  "6",  NULL, NULL, },
        { NULL, "1",  "2",  "3",  NULL, NULL, },
        { NULL, NULL, "0",  ".",  NULL, NULL, },
        { NULL, NULL, NULL, NULL, NULL, NULL, },
    };

    static const char *macpad[KEYPAD_Y_SIZE][KEYPAD_X_SIZE] = {
        { NULL, NULL, NULL, NULL, NULL, NULL, },
        { NULL, NULL, "=",  "/",  "*",  NULL, },
        { NULL, "7",  "8",  "9",  "-",  NULL, },
        { NULL, "4",  "5",  "6",  "+",  NULL, },
        { NULL, "1",  "2",  "3",  NULL, NULL, },
        { NULL, "0",  ".",  NULL, NULL, NULL, },
        { NULL, NULL, NULL, NULL, NULL, NULL, },
    };

    static const char *alphabet[ALPHABET_Y_SIZE][ALPHABET_X_SIZE] = {
        { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
        { NULL, "aA", "bB", "cC", "dD", "eE", "fF", "gG", "hH", "iI", "jJ", "kK", "lL", "mM", "nN", "oO", "pP", "qQ", "rR", "sS", "tT", "uU", "vV", "wW", "xX", "yY", "zZ", NULL },
        { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
    };

    make_spatial_graph_iter(&zxcvbn->spatial_graph_qwerty, (const char **)qwerty,
                            KEYBRD_X_SIZE, KEYBRD_Y_SIZE, 2, 6, get_slant_coords);

    make_spatial_graph_iter(&zxcvbn->spatial_graph_dvorak, (const char **)dvorak,
                            KEYBRD_X_SIZE, KEYBRD_Y_SIZE, 2, 6, get_slant_coords);

    make_spatial_graph_iter(&zxcvbn->spatial_graph_keypad, (const char **)keypad,
                            KEYPAD_X_SIZE, KEYPAD_Y_SIZE, 1, 8, get_align_coords);

    make_spatial_graph_iter(&zxcvbn->spatial_graph_macpad, (const char **)macpad,
                            KEYPAD_X_SIZE, KEYPAD_Y_SIZE, 1, 8, get_align_coords);

    make_spatial_graph_iter(&zxcvbn->spatial_graph_alphabet, (const char **)alphabet,
                            ALPHABET_X_SIZE, ALPHABET_Y_SIZE, 2, 2, get_liner_coords);

#undef KEYBRD_X_SIZE
#undef KEYBRD_Y_SIZE
#undef KEYPAD_X_SIZE
#undef KEYPAD_Y_SIZE
#undef ALPHABET_X_SIZE
#undef ALPHABET_Y_SIZE
}

void
zxcvbn_res_init(struct zxcvbn_res *res, struct zxcvbn *zxcvbn)
{
    res->zxcvbn = zxcvbn;
    res->n_matches = 0;
    res->matches = res->match_buf;
    res->n_matches_reserved = ARRAY_SIZE(res->match_buf);
}

void
zxcvbn_res_release(struct zxcvbn_res *res)
{
    if (res->matches != res->match_buf)
        __free(res->zxcvbn, res->matches);
}

static struct zxcvbn_match *
match_add(struct zxcvbn_res *res)
{
    size_t size;
    if (res->n_matches_reserved == res->n_matches) {
        res->n_matches_reserved += ARRAY_SIZE(res->match_buf);
        size = sizeof(struct zxcvbn_match) * res->n_matches_reserved;
        if (res->matches == res->match_buf) {
            if ((res->matches = __malloc(res->zxcvbn, size)) == NULL)
                return NULL;
            memcpy(res->matches, res->match_buf, sizeof(res->match_buf));
        } else {
            if ((res->matches = __realloc(res->zxcvbn, res->matches, size)) == NULL)
                return NULL;
        }
    }

    return res->matches + res->n_matches++;
}

static struct zxcvbn_match *
push_match(struct zxcvbn_res *res,
           enum zxcvbn_match_type type, struct zxcvbn_spatial_graph *spatial_graph,
           unsigned int i, unsigned int j, unsigned int turns, unsigned int shifted)
{
    struct zxcvbn_match *match;

    if ((match = match_add(res)) == NULL)
        return NULL;

    match->type = type;
    match->spatial_graph = spatial_graph;
    match->i = i;
    match->j = j;
    match->turns = turns;
    match->shifted = shifted;

    return match;
}

static struct zxcvbn_match *
push_match_date(struct zxcvbn_res *res,
                unsigned int i, unsigned int j, unsigned int separator, unsigned int long_year, unsigned int only_year)
{
    struct zxcvbn_match *match;

    if ((match = match_add(res)) == NULL)
        return NULL;
    
    match->type = ZXCVBN_MATCH_TYPE_DATE;
    match->i = i;
    match->j = j;
    match->separator = separator;
    match->long_year = long_year;
    
    return match;
}

static struct zxcvbn_match *
push_match_dict(struct zxcvbn_res *res, unsigned int i, unsigned int j, unsigned int rank)
{
    struct zxcvbn_match *match;

    if ((match = match_add(res)) == NULL)
        return NULL;

    match->type = ZXCVBN_MATCH_TYPE_DICT;
    match->i = i;
    match->j = j;
    match->rank = rank;

    return match;
}

static struct zxcvbn_match *
push_match_bruteforce(struct zxcvbn_res *res,
                      unsigned int i, unsigned int j, unsigned int bruteforce_card)
{
    struct zxcvbn_match *match;

    if ((match = match_add(res)) == NULL)
        return NULL;

    match->type = ZXCVBN_MATCH_TYPE_BRUTEFORCE;
    match->i = i;
    match->j = j;
    match->entropy = log2(pow(bruteforce_card, j - i + 1));

    return match;
}

int
match_spatial_iter(struct zxcvbn_res *res,
                   const char *password, unsigned int password_len, struct zxcvbn_spatial_graph *spatial_graph)
{
    int i, j, k, cur_dir, prv_dir, turns, shifted;
    char cur, prv;
    const char *s, *p;

    i = j = 0;
    prv_dir = -1;
    turns = 0;
    shifted = 0;
    
    while (1) {
        prv = password[j];
        ++j;
        cur = password[j];

        for (cur_dir = 0; cur_dir < spatial_graph->n_coords; ++cur_dir) {
            s = spatial_graph->data[prv][cur_dir];
            for (p = s; p - s < spatial_graph->token_size; ++p) {
                if (*p == cur) {
                    shifted += p != s;
                    goto found;
                }
            }
        }

        if (j - i > 2) {
            if (push_match(res, ZXCVBN_MATCH_TYPE_SPATIAL, spatial_graph, i, j - 1, turns, shifted) == NULL)
                return -1;
        }

        i = j;
        prv_dir = -1;
        turns = 0;
        shifted = 0;

        if (i == password_len)
            break;

        continue;
found:
        if (cur_dir != prv_dir) {
            ++turns;
            prv_dir = cur_dir;
        }
    }
 
    return 0;
}

static int
match_spatial(struct zxcvbn_res *res, const char *password, unsigned int password_len)
{
    struct zxcvbn *zxcvbn;

    zxcvbn = res->zxcvbn;

    if (match_spatial_iter(res, password, password_len, &zxcvbn->spatial_graph_qwerty) ||
        match_spatial_iter(res, password, password_len, &zxcvbn->spatial_graph_dvorak) ||
        match_spatial_iter(res, password, password_len, &zxcvbn->spatial_graph_keypad) ||
        match_spatial_iter(res, password, password_len, &zxcvbn->spatial_graph_macpad) ||
        match_spatial_iter(res, password, password_len, &zxcvbn->spatial_graph_alphabet))
        return -1;

    return 0;
}

static int
match_digits(struct zxcvbn_res *res, const char *password, unsigned int password_len)
{
    int i, j;

    i = j = 0;

    while (1) {
        if (isdigit(password[j])) {
            ++j;
        } else {
            if (j - i > 2) {
                if (push_match(res, ZXCVBN_MATCH_TYPE_DIGITS, NULL, i, j - 1, 0, 0) == NULL)
                    return -1;
            }
            if (j == password_len)
                break;
            ++j;
            i = j;
        }
    }

    return 0;
}

static inline int
date_is_yyyy(const char *yyyy)
{
    switch (*yyyy) {
    case '1':
        return *(yyyy + 1) == '9';
    case '2':
        return *(yyyy + 1) == '0' && (*(yyyy + 2) == '0' || *(yyyy + 2) == '1');
    }

    return 0;
}

static inline int
date_is_yy(const char *yy)
{
    return isdigit(*yy) && isdigit(*(yy + 1));
}

static inline int
date_is_dd(const char *dd)
{
    unsigned int dd_10;

    dd_10 = 10*(*dd - '0') + *(dd + 1) - '0';

    return dd_10 <= 31;
}

static inline int
date_is_mm(const char *mm)
{
    unsigned int mm_10;

    mm_10 = 10*(*mm - '0') + *(mm + 1) - '0';

    return mm_10 <= 12;
}

static int
date_nosep(const char **date, const char *buf, unsigned int len)
{
    int i, is_yyyy, date_len, is_ddmm, is_mmdd;

    date_len = -1;

    for (i = 0; i < len - 3; ++i) {

        is_yyyy = date_is_yyyy(buf + i);
        is_ddmm = date_is_dd(buf + i) && date_is_mm(buf + i + 2);
        is_mmdd = date_is_mm(buf + i) && date_is_dd(buf + i + 2);

        if (len - i >= 8) {
            if (is_yyyy) {
                if ((date_is_mm(buf + i + 4) && date_is_dd(buf + i + 4 + 2)) ||
                    (date_is_dd(buf + i + 4) || date_is_mm(buf + i + 4 + 2))) {
                    date_len = 8;
                    break;
                }
            }

            if (date_is_yyyy(buf + i + 4)) {
                if (is_mmdd || is_ddmm) {
                    date_len = 8;
                    break;
                }
            }
        }

        if (len - i >= 6) {
            if (date_is_yy(buf + i)) {
                if ((date_is_mm(buf + i + 2) && date_is_dd(buf + i + 2 + 2)) ||
                    (date_is_dd(buf + i + 2) && date_is_mm(buf + i + 2 + 2))) {
                    date_len = 6;
                    break;
                }
            }

            if (date_is_yy(buf + i + 4)) {
                 if (is_mmdd || is_ddmm) {
                    date_len = 6;
                    break;
                }
            }
        }

        if (is_yyyy) {
            date_len = 4;
            break;
        }
    }

    if (date_len > 0)
        *date = buf + i;

    return date_len;
}

static int
match_date(struct zxcvbn_res *res, const char *password, int password_len)
{
    int i, j, maybe_full_year, maybe_separator, date_len, long_year, only_year;
    unsigned int len;
    const char *date;

    for (i = 0; i < password_len - 5; i += len) {
        len = 0;

        if (!isdigit(password[i + len++]) || !isdigit(password[i + len++]))
            continue;

        for (j = i + len; isdigit(password[j]); ++j)
            ++len;

        maybe_full_year = password_len - i >= 10;
        maybe_separator = password_len - i >= 8;

        switch (len) {
        case 1:
        case 3:
        case 5:
            break;

        case 2:
            if (maybe_separator == 0 ||
                !((date_is_mm(password + i) && date_is_dd(password + i + 2 + 1)) ||
                  (date_is_dd(password + i) && date_is_mm(password + i + 2 + 1)))) {
                break;
            }

            if (maybe_full_year == 1) {
                if (date_is_yyyy(password + i + 2 + 1 + 2 + 1)) {
                    len = 10;
                    if (push_match_date(res, i, i + len - 1, 1, 1, 0) == NULL)
                        return -1;
                }
            } else {
                if (date_is_yy(password + i + 2 + 1 + 2 + 1)) {
                    len = 8;
                    if (push_match_date(res, i, i + len - 1, 1, 0, 0) == NULL)
                        return -1;
                }
            }

            break;

        case 4:
            if (maybe_full_year == 0 || !date_is_yyyy(password + i))
                break;            

            if ((date_is_mm(password + i + 4 + 1) && date_is_dd(password + i + 4 + 1 + 2 + 1)) ||
                (date_is_dd(password + i + 4 + 1) && date_is_mm(password + i + 4 + 1 + 2 + 1))) {
                len = 10;
                only_year = 0;
            } else {
                only_year = 1;
            }

            if (push_match_date(res, i, i + len - 1, 1, 1, only_year) == NULL)
                return -1;

            break;

        default:
            if ((date_len = date_nosep(&date, password + i, len)) > 0) {
                long_year = 1;
                only_year = 0;
                switch (date_len) {
                case 4:
                    only_year = 1;
                    break;
                case 6:
                    long_year = 0;
                    break;
                case 8:
                    break;
                default:
                    assert(0);
                }

                i += date - (password + i);
                len = date_len;

                if (push_match_date(res, i, i + len - 1, 0, long_year, only_year) == NULL)
                    return -1;
            }

            break;
        }
    }

    return 0;
}

static int
match_dict_iter(struct zxcvbn_res *res, struct zxcvbn_dict *dict, const char *password, unsigned int password_len)
{
    int i, j;
    struct zxcvbn_node *node, *parent;

    for (i = 0; i < password_len; ++i) {
        parent = dict->root;
        for (j = i; j < password_len; ++j) {
            if ((node = parent->children[password[j]]) == NULL)
                break;
            if (node->rank > 0) {
                if (push_match_dict(res, i, j, node->rank) == NULL)
                    return -1;
            }
            parent = node;
        }
    }

    return 0;
}

static char *
pack_word(struct zxcvbn *zxcvbn, char *dst, const char *src, unsigned int len)
{
    int i;

    for (i = 0; i < len; ++i)
        dst[i] = zxcvbn->pack_table[tolower(src[i])];

    return dst;
}

static int
match_dict(struct zxcvbn_res *res, const char *password, unsigned int password_len,
           char **dict_words, unsigned int n_dict_words)
{
    int i, dict_word_len;
    char pack_password[ZXCVBN_PASSWORD_LEN_MAX], *pack_dict_word, *s;
    struct zxcvbn_dict *dict;
    struct zxcvbn *zxcvbn;

    zxcvbn = res->zxcvbn;

    pack_word(zxcvbn, pack_password, password, password_len);

    for (i = 0; i < n_dict_words; ++i) {
        dict_word_len = strlen(dict_words[i]);
        pack_dict_word = pack_word(zxcvbn, dict_words[i], dict_words[i], dict_word_len);

        if ((s = memmem(pack_password, password_len, pack_dict_word, dict_word_len)) != NULL) {
            push_match_dict(res, s - pack_password, s - pack_password + dict_word_len - 1, 1);
        }
    }

    LIST_FOREACH(dict, &res->zxcvbn->dict_head, list) {
        if (match_dict_iter(res, dict, pack_password, password_len) < 0)
            return -1;
    }

    return 0;
}

static unsigned int
nCk(unsigned int n, unsigned int k)
{
    unsigned int r, d;

    if (k > n)
        return 0;
    if (k == 0)
        return 1;

    r = 1;

    for (d = 1; d <= k; ++d) {
        r *= n;
        r /= d;
        n -= 1;
    }

    return r;
}

static void
entropy_dict(struct zxcvbn *zxcvbn, struct zxcvbn_match *match, const char *password, unsigned int password_len)
{
    int i, ch, upper, lower, min_lower_upper;
    unsigned long possibilities;

    match->entropy = log2(match->rank);

    upper = 0;
    lower = 0;

    for (i = match->i; i <= match->j; ++i) {
        ch = password[i];
        if (isalpha(ch)) {
            if (isupper(ch))
                ++upper;
            else if (islower(ch))
                ++lower;
        }
    }

    if (upper) {
        printf("XXXX %d\n", upper);
        min_lower_upper = MIN(lower, upper);
        possibilities = 0;
        for (i = 0; i <= min_lower_upper; ++i)
            possibilities += nCk(upper + lower, i);
        match->entropy += log2(possibilities);
    }
}

static void
entropy_spatial(struct zxcvbn *zxcvbn, struct zxcvbn_match *match)
{
    int i, j;
    unsigned int length, turns, possible_turns, S, U, min_SU;
    unsigned long possibilities;
    struct zxcvbn_spatial_graph *spatial_graph;

    spatial_graph = match->spatial_graph;

    length = match->j - match->i + 1;
    turns = match->turns;
    possibilities = 0;

    for (i = 2; i <= length; ++i) {
        possible_turns = MIN(turns, i - 1);
        for (j = 1; j <= possible_turns; ++j) {
            possibilities += nCk(i - 1, j - 1) * spatial_graph->n_chars * pow(spatial_graph->degree, j);
        }
    }

    match->entropy = log2(possibilities);

    if (match->shifted) {
        S = match->shifted;
        U = length - S;
        min_SU = MIN(S, U);
        possibilities = 0;
        for (i = 0; i <= min_SU; ++i)
            possibilities += nCk(S + U, i);
        match->entropy += log2(possibilities);
    }
}

static void
entropy_digits(struct zxcvbn *zxcvbn, struct zxcvbn_match *match)
{
    match->entropy = log2(pow(10, match->j - match->i + 1));
}

static void
entropy_date(struct zxcvbn *zxcvbn, struct zxcvbn_match *match)
{
    unsigned int num_years;
    unsigned int num_days;

    num_years = match->long_year ? 119 : 100;
    num_days = match->only_year ? 1 : 31 * 12;
    match->entropy = log2(num_days * num_years);

    if (match->separator)
        match->entropy += 2;
}

static int
calc_bruteforce_card(const char *password, unsigned int password_len, unsigned int n_symbols)
{
    int i, digit, lower, upper, symbol;

    digit = lower = upper = symbol = 0;

    for (i = 0; i < password_len; ++i) {
        switch (password[i]) {
        case '0'...'9':
            digit = 1;
            break;
        case 'a'...'z':
            lower = 1;
            break;
        case 'A'...'Z':
            upper = 1;
            break;
        default:
            symbol = 1;
            break;
        }
    }

    return digit * ('9' - '0') +
           lower * ('z' - 'a') +
           upper * ('Z' - 'A') + symbol * n_symbols;
}



struct zxcvbn *
zxcvbn_init(struct zxcvbn *zxcvbn_buf,
            void *(*zxcvbn_malloc)(size_t size),
            void *(*zxcvbn_realloc)(void *ptr, size_t size),
            void (*zxcvbn_free)(void *ptr),
            const char *symbols)
{
    int i, l33t, ch;
    const char *s;
    struct zxcvbn *zxcvbn;

    if (zxcvbn_malloc == NULL) {
        assert(zxcvbn_realloc == NULL);
        assert(zxcvbn_free == NULL);

        zxcvbn_malloc = malloc;
        zxcvbn_realloc = realloc;
        zxcvbn_free = free;
    } else {
        assert(zxcvbn_realloc != NULL);
        assert(zxcvbn_free != NULL);
    }

    if (zxcvbn_buf == NULL) {
        if ((zxcvbn = (*zxcvbn_malloc)(sizeof(struct zxcvbn))) == NULL)
            return NULL;
    } else {
        zxcvbn = zxcvbn_buf;
    }

    memset(zxcvbn, 0, sizeof(*zxcvbn));

    LIST_INIT(&zxcvbn->dict_head);

    zxcvbn->zxcvbn_malloc = zxcvbn_malloc;
    zxcvbn->zxcvbn_realloc = zxcvbn_realloc;
    zxcvbn->zxcvbn_free = zxcvbn_free;

    zxcvbn->allocated = zxcvbn_buf == NULL;

    memset(zxcvbn->pack_table, '.', sizeof(zxcvbn->pack_table));

    for (i = 'a'; i <= 'z'; ++i)
        zxcvbn->pack_table[i] = zxcvbn->pack_table_size++;

    for (i = '0'; i <= '9'; ++i)
        zxcvbn->pack_table[i] = zxcvbn->pack_table_size++;

    for (s = symbols; *s != '\0'; ++s) {
        if (zxcvbn->pack_table[*s] == '.') {
            zxcvbn->pack_table[*s] = zxcvbn->pack_table_size++;
            zxcvbn->n_symbols++;
        }
    }

    // l33t
    for (i = 0; i < 256; ++i) { 
        ch = zxcvbn->pack_table[i];
        switch (ch) {
        case '4':
        case '@':
            l33t = 'a';
            break;
        case '$':
        case '5':
            l33t = 's';
            break;
        default:
            continue;
        }
        zxcvbn->pack_table_size--;
        zxcvbn->pack_table[ch] = zxcvbn->pack_table[l33t];
    }

    make_spatial_graph(zxcvbn);

    return zxcvbn;
}

static int
min_entropy(struct zxcvbn_res *res, const char *password, unsigned int password_len)
{
    int pos, match_i;
    double pos_entropy[ZXCVBN_PASSWORD_LEN_MAX], entropy;
    unsigned int bruteforce_card;
    struct zxcvbn_match *matches[ZXCVBN_PASSWORD_LEN_MAX], *match, *match_bruteforce;
    struct zxcvbn_match_head *match_head;

    bruteforce_card = calc_bruteforce_card(password, password_len, res->zxcvbn->n_symbols);
    pos_entropy[0] = 0;

    for (pos = 0; pos < password_len; ++pos) {
        pos_entropy[pos] = pos > 0 ? pos_entropy[pos - 1] : 0;
        pos_entropy[pos] += log2(bruteforce_card);
        matches[pos] = NULL;

        for (match_i = 0; match_i < res->n_matches; ++match_i) {
            match = res->matches + match_i;
            if (match->j != pos)
                continue;

            entropy = match->i > 0 ? pos_entropy[match->i - 1] : 0;
            entropy += match->entropy;

            if (pos_entropy[pos] > entropy) {
                pos_entropy[pos] = entropy;
                matches[pos] = match;
            }
        }
    }

    res->entropy = pos_entropy[password_len - 1];

    match_head = &res->match_head;
    CIRCLEQ_INIT(match_head);

    pos = password_len - 1;
    while (pos >= 0) {
        if ((match = matches[pos]) == NULL) {
            --pos;
        } else {
            CIRCLEQ_INSERT_HEAD(match_head, match, list);
            pos = match->i - 1;
        }
    }

    pos = 0;
    CIRCLEQ_FOREACH(match, match_head, list) {
        if (match->i > pos) {
            if ((match_bruteforce = push_match_bruteforce(res, pos, match->i - 1, bruteforce_card)) == NULL)
                return -1;
            CIRCLEQ_INSERT_BEFORE(match_head, match, match_bruteforce, list);
        }
        pos = match->j + 1;
    }

    if (pos < password_len) {
        if ((match_bruteforce = push_match_bruteforce(res, pos, password_len - 1, bruteforce_card)) == NULL)
            return -1;
        CIRCLEQ_INSERT_TAIL(match_head, match_bruteforce, list);
    }

    return 0;
}

int
zxcvbn_match(struct zxcvbn_res *res, const char *password, unsigned int password_len,
             char **dict_words, unsigned int n_dict_words)
{
    int i;
    char pack_password[ZXCVBN_PASSWORD_LEN_MAX];
    struct zxcvbn *zxcvbn;
    struct zxcvbn_match *match;

    assert(password_len > 0);
    assert(password_len <= ZXCVBN_PASSWORD_LEN_MAX);

    if (match_spatial(res, password, password_len))
        return -1;
    if (match_digits(res, password, password_len))
        return -1;
    if (match_date(res, password, password_len))
        return -1;
    if (match_dict(res, password, password_len, dict_words, n_dict_words))
        return -1;

    zxcvbn = res->zxcvbn;

    for (i = 0; i < res->n_matches; ++i) {
        match = res->matches + i;
        switch (match->type) {
        case ZXCVBN_MATCH_TYPE_DICT:
            entropy_dict(zxcvbn, match, password, password_len);
            break;        
        case ZXCVBN_MATCH_TYPE_SPATIAL:
            entropy_spatial(zxcvbn, match);
            break;
        case ZXCVBN_MATCH_TYPE_DIGITS:
            entropy_digits(zxcvbn, match);
            break;
        case ZXCVBN_MATCH_TYPE_DATE:
            entropy_date(zxcvbn, match);
            break;
        default:
            assert(0);
        }
    }

    return min_entropy(res, password, password_len);
}

static void
free_node(struct zxcvbn *zxcvbn, struct zxcvbn_node *node)
{
    int i;
    struct zxcvbn_node *child;

    for (i = 0; i < zxcvbn->pack_table_size; ++i) {
        if ((child = node->children[i]) != NULL)
            free_node(zxcvbn, child);
    }

    __free(zxcvbn, node);
}

void
zxcvbn_dict_release(struct zxcvbn_dict *dict)
{
    LIST_REMOVE(dict, list);

    free_node(dict->zxcvbn, dict->root);

    if (dict->allocated)
        __free(dict->zxcvbn, dict);
}

void
zxcvbn_release(struct zxcvbn *zxcvbn)
{
    struct zxcvbn_dict *dict;

    if (zxcvbn == NULL)
        return;

    while (!LIST_EMPTY(&zxcvbn->dict_head)) {
        dict = LIST_FIRST(&zxcvbn->dict_head);
        zxcvbn_dict_release(dict);
    }

    if (zxcvbn->allocated)
        __free(zxcvbn, zxcvbn);
}

static struct zxcvbn_node *
make_node(struct zxcvbn *zxcvbn)
{
    unsigned int size;
    char *memory;
    struct zxcvbn_node *node;

    size = sizeof(*node) + zxcvbn->pack_table_size * sizeof(struct zxcvbn_node *);

    if ((memory =__malloc(zxcvbn, size)) == NULL)
        return NULL;

    memset(memory, 0, size);
    node = (struct zxcvbn_node *)memory;
    node->children = (struct zxcvbn_node **)(memory + sizeof(*node));
    node->rank = -1;

    return node;
}

struct zxcvbn_dict *
zxcvbn_dict_init(struct zxcvbn *zxcvbn, struct zxcvbn_dict *dict_buf, const char *name)
{
    struct zxcvbn_dict *dict;

    if (dict_buf != NULL) {
        dict = dict_buf;
        dict->allocated = 0;
    } else {
        if ((dict = __malloc(zxcvbn, sizeof(struct zxcvbn_dict))) == NULL)
            return NULL;
        dict->allocated = 1;
    }
        
    dict->zxcvbn = zxcvbn;
    strncpy(dict->name, name, sizeof(dict->name));
    dict->name[sizeof(dict->name) - 1];
    if ((dict->root = make_node(zxcvbn)) == NULL) {
        zxcvbn_dict_release(dict);
        return NULL;
    }

    LIST_INSERT_HEAD(&zxcvbn->dict_head, dict, list);

    return dict;
}

int
zxcvbn_dict_add_word(struct zxcvbn_dict *dict, const char *word, unsigned int word_len, unsigned int rank)
{
    int i;
    unsigned int len;
    const char *s;
    const void *ptr;
    char word_buf[ZXCVBN_PASSWORD_LEN_MAX];
    struct zxcvbn_node *node, *parent;

    len = MIN(word_len, sizeof(word_buf));

    if (pow(26, len) < rank) {
        // bruteforce possibilities are less then word rank.
        return 0;
    }

    pack_word(dict->zxcvbn, word_buf, word, len);

    parent = dict->root;

    for (i = 0;; ++i) {
        if ((node = parent->children[word_buf[i]]) == NULL) {
            if ((node = make_node(dict->zxcvbn)) == NULL)
                return -1;
            parent->children[word_buf[i]] = node; 
        }

        if (i == len - 1) {
            if (node->rank == -1)
                node->rank = rank;
            break;
        }

        parent = node;
    }

    return 1;
}
