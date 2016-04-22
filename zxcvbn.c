#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

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
    return digit * ('9' - '0' + 1) +
           lower * ('z' - 'a' + 1) +
           upper * ('Z' - 'A' + 1) + symbol * n_symbols;
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
    case ZXCVBN_MATCH_TYPE_SEQUENCE:
        return "sequence";
    case ZXCVBN_MATCH_TYPE_REPEAT:
        return "repeat";
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

static void
make_spatial_graph(struct zxcvbn *zxcvbn)
{
#define KEYBRD_X_SIZE 16
#define KEYBRD_Y_SIZE 6
#define KEYPAD_X_SIZE 6
#define KEYPAD_Y_SIZE 7

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

    make_spatial_graph_iter(&zxcvbn->spatial_graph_qwerty, (const char **)qwerty,
                            KEYBRD_X_SIZE, KEYBRD_Y_SIZE, 2, 6, get_slant_coords);

    make_spatial_graph_iter(&zxcvbn->spatial_graph_dvorak, (const char **)dvorak,
                            KEYBRD_X_SIZE, KEYBRD_Y_SIZE, 2, 6, get_slant_coords);

    make_spatial_graph_iter(&zxcvbn->spatial_graph_keypad, (const char **)keypad,
                            KEYPAD_X_SIZE, KEYPAD_Y_SIZE, 1, 8, get_align_coords);

    make_spatial_graph_iter(&zxcvbn->spatial_graph_macpad, (const char **)macpad,
                            KEYPAD_X_SIZE, KEYPAD_Y_SIZE, 1, 8, get_align_coords);

#undef KEYBRD_X_SIZE
#undef KEYBRD_Y_SIZE
#undef KEYPAD_X_SIZE
#undef KEYPAD_Y_SIZE
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
        if (!res->zxcvbn->max_matches_num)
            res->n_matches_reserved += ARRAY_SIZE(res->match_buf);
        else {
            if (res->n_matches_reserved >= res->zxcvbn->max_matches_num)
                return NULL;
            res->n_matches_reserved =
                    MIN(res->n_matches_reserved + ARRAY_SIZE(res->match_buf),
                        res->zxcvbn->max_matches_num);
        }
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

static int
match_spatial_iter(struct zxcvbn_res *res,
                   const char *password, unsigned int password_len, struct zxcvbn_spatial_graph *spatial_graph)
{
    int i, j, cur_dir, prv_dir, turns, shifted;
    unsigned char cur, prv;
    const char *s, *p;

    i = j = 0;
    prv_dir = -1;
    turns = 0;
    shifted = 0;
    
    while (i + 2 < password_len) {
        prv = password[j];
        ++j;
        if (j < password_len) {
            cur = password[j];

            for (cur_dir = 0; cur_dir < spatial_graph->n_coords; ++cur_dir) {
                s = spatial_graph->data[prv][cur_dir];
                for (p = s; p - s < spatial_graph->token_size; ++p) {
                    if (*p == cur) {
                        shifted += p != s;
                        if (cur_dir != prv_dir) {
                            ++turns;
                            prv_dir = cur_dir;
                        }
                        continue;
                    }
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
        match_spatial_iter(res, password, password_len, &zxcvbn->spatial_graph_macpad))
        return -1;

    return 0;
}

/* Repeat =================================================================== */

static int8_t
zxcvbn_repeat_match(struct zxcvbn_res *res,
                    char *password, uint32_t password_len)
{
    char ch;
    uint32_t i, j;

    i = 0;
    while (i + 1 < password_len) {
        ch = password[i];
        for (j = i + 1; j < password_len && password[j] == ch; j++)
            {;}
        if (j - i > 2) {
            if (!push_match(res, ZXCVBN_MATCH_TYPE_REPEAT,
                            NULL, i, j - 1, 0, 0))
                return -1;
        }
        i = j;
    }
    return 0;
}

static void
zxcvbn_repeat_calculate_entropy(struct zxcvbn *zxcvbn,
                                struct zxcvbn_match *match,
                                char *password, uint32_t password_len)
{
    match->entropy = log2(calc_bruteforce_card(password + match->i, 1,
                                               zxcvbn->n_symbols) *
                          (match->j - match->i + 1));
}

/* Repeat ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

/* Sequence ================================================================= */

#define ZXCVBN_SEQUENCES_DEF(m) \
    m("abcdefghijklmnopqrstuvwxyz", 0)          \
    m("ABCDEFGHIJKLMNOPQRSTUVWXYZ", 1)          \
    m("f,dult;pbqrkvyjghcnea[wxio]sm'.z", 1)    \
    m("F<DULT:PBQRKVYJGHCNEA{WXIO}SM\">Z", 2)   \
    m("abvgdegziyklmnoprstufhc", 1)             \
    m("ABVGDEGZIYKLMNOPRSTUFHC", 2)             \
    m("0123456789", 0)                          \

#define ZXCVBN_SEQUENCE_OBVIOUS_START   "aAzZfF019"
#define ZXCVBN_SEQUENCE_MIN_LEN         3

struct zxcvbn_sequence {
    char           *str;
    unsigned int    len;
    unsigned int    extra_entropy;
};

static struct zxcvbn_sequence zxcvbn_sequences[] = {
#define ZXCVBN_M(s,e)   {s, sizeof(s) - 1, e},
    ZXCVBN_SEQUENCES_DEF(ZXCVBN_M)
#undef ZXCVBN_M
};

static struct zxcvbn_match *
zxcvbn_sequence_add_match(struct zxcvbn_res *res, uint32_t i, uint32_t j,
                          struct zxcvbn_sequence *seq, int8_t dir)
{
    struct zxcvbn_match *match;

    if (!(match = match_add(res)))
        return NULL;

    match->type = ZXCVBN_MATCH_TYPE_SEQUENCE;
    match->i = i;
    match->j = j;
    match->seq = seq;
    match->flags |= (dir == -1 ? ZXCVBN_MATCH_DESC_SEQ : 0);
    return match;
}

static int8_t
zxcvbn_sequence_match(struct zxcvbn_res *res,
                      char *password, uint32_t password_len)
{
    struct zxcvbn_sequence *seq;
    uint32_t i, j, s, i_n, j_n, prev_n;
    int32_t dir;
    char *p;

    i = 0;
    while (i + ZXCVBN_SEQUENCE_MIN_LEN - 1 < password_len) {
        j = i + 1;
        for (s = 0; s < ARRAY_SIZE(zxcvbn_sequences); s++) {
            seq = zxcvbn_sequences + s;
            if (!(p = strchr(seq->str, password[i])))
                continue;
            i_n = p - seq->str;
            if (!(p = strchr(seq->str, password[j])))
                continue;
            j_n = p - seq->str;
            if ((i_n + 1) % seq->len == j_n) {
                dir = 1;
                break;
            } else if ((j_n + 1) % seq->len == i_n) {
                dir = -1;
                break;
            }
        }
        if (s == ARRAY_SIZE(zxcvbn_sequences)) {
            i++;
            continue;
        }

        j++;
        while (j < password_len) {
            p = strchr(seq->str, password[j]);
            if (!p)
                break;
            prev_n = j_n;
            j_n = p - seq->str;
            if (j_n != (seq->len + prev_n + dir) % seq->len)
                break;
            j++;
        }

        if (j - i >= ZXCVBN_SEQUENCE_MIN_LEN) {
            if (!zxcvbn_sequence_add_match(res, i, j, seq, dir))
                return -1;
        }
        i = j;
    }
    return 0;
}

static void
zxcvbn_sequence_calculate_entropy(struct zxcvbn *zxcvbn,
                                  struct zxcvbn_match *match,
                                  char *password, uint32_t password_len)
{
    if (strchr(ZXCVBN_SEQUENCE_OBVIOUS_START, password[match->i]))
        match->entropy = 1;
    else
        match->entropy = log2(match->seq->len) + match->seq->extra_entropy;
    if (match->flags & ZXCVBN_MATCH_DESC_SEQ)
        match->entropy++;
    match->entropy += log2(match->j - match->i + 1);
}

/* Sequence ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

static int
match_digits(struct zxcvbn_res *res,
             const char *password, unsigned int password_len)
{
    int last = -1, i;

    for (i = 0; i <= password_len; i++) {
        if (i == password_len || !isdigit(password[i])) {
            if (i - last - 1 > 2 &&
                    !push_match(res, ZXCVBN_MATCH_TYPE_DIGITS, NULL,
                                last + 1, i - 1, 0, 0))
                return -1;
            last = i;
        }
    }
    return 0;
}

/* Date ===================================================================== */

#define ZXCVBN_DATE_MIN_NOSEP_LEN       4
#define ZXCVBN_DATE_MAX_NOSEP_LEN       8
#define ZXCVBN_DATE_MIN_SEP_LEN         6
#define ZXCVBN_DATE_MAX_SEP_LEN         10

#define ZXCVBN_DATE_REF_YEAR            2000
#define ZXCVBN_DATE_MIN_YEAR            1000
#define ZXCVBN_DATE_MAX_YEAR            2050
#define ZXCVBN_DATE_MIN_YEAR_SPACE      20

#define ZXCVBN_DATE_PROBE_LEFT_YEAR     (1 << 0)
#define ZXCVBN_DATE_PROBE_RIGHT_YEAR    (1 << 1)
#define ZXCVBN_DATE_PROBE_FULL_YEAR     (1 << 2)

struct zxcvbn_date_state {
    int8_t      next[3];
    uint8_t     skip[3];
    int8_t      num;
    uint8_t     try;
    uint8_t     probe_flags;
};

static uint16_t
zxcvbn_parse_number(char *str, uint8_t len)
{
    uint16_t n = 0;

    while (len--)
        n = n * 10 + *str++ - '0';
    return n;
}

static inline uint8_t
zxcvbn_date_probe_year(struct zxcvbn_date *date, char *str)
{
    uint16_t year;

    year = zxcvbn_parse_number(str, 4);
    if (year < ZXCVBN_DATE_MIN_YEAR ||
            year > ZXCVBN_DATE_MAX_YEAR)
        return 0;
    memset(date, 0, sizeof(*date));
    date->year = year;
    date->flags = ZXCVBN_DATE_ONLY_YEAR | ZXCVBN_DATE_FULL_YEAR;
    return 1;
}

static uint8_t
zxcvbn_date_probe(struct zxcvbn_date *date, uint16_t *nums, uint8_t flags,
                  struct zxcvbn_date *dates, uint32_t dates_num)
{
    static uint8_t meanings[] = {2, 1, 1, 2, 0, 1, 1, 0};
    struct zxcvbn_date temp, best;
    uint8_t over_31 = 0, over_12 = 0, equal_0 = 0;
    uint8_t m, i, j, *p;

    if (nums[1] > 31 || nums[1] == 0)
        return 0;
    for (i = 0; i < 3; i++) {
        if (nums[i] > 31)
            over_31++;
        if (nums[i] > 12)
            over_12++;
        if (nums[i] == 0)
            equal_0++;
    }
    if (over_31 >= 2 || over_12 == 3 || equal_0 >= 2)
        return 0;

    memset(&best, 0, sizeof(best));
    for (m = 1; m < 4; m <<= 1) {
        if (!(flags & m))
            continue;
        memset(&temp, 0, sizeof(temp));
        temp.flags = (flags & ZXCVBN_DATE_PROBE_FULL_YEAR ?
                      ZXCVBN_DATE_FULL_YEAR : 0);
        temp.year = nums[m & 2];
        if (temp.flags & ZXCVBN_DATE_FULL_YEAR) {
            if ((temp.year < ZXCVBN_DATE_MIN_YEAR ||
                    temp.year > ZXCVBN_DATE_MAX_YEAR))
                continue;
        } else
            temp.year += ZXCVBN_DATE_REF_YEAR - (temp.year > 50 ? 100 : 0);
        p = meanings + ((m & 2) << 1);
        for (j = 0; j < 2; j++, p += 2) {
            temp.day   = nums[p[0]];
            temp.month = nums[p[1]];
            if (temp.day == 0 || temp.day > 31 ||
                    temp.month == 0 || temp.month > 12)
                continue;
            for (i = 0; i < dates_num; i++) {
                if (temp.day == dates[i].day &&
                        temp.month == dates[i].month &&
                        temp.year == dates[i].year) {
                    memcpy(date, &temp, sizeof(temp));
                    date->flags |= ZXCVBN_DATE_FROM_LIST;
                    return 1;
                }
            }
            if (!best.day ||
                    abs(best.year - ZXCVBN_DATE_REF_YEAR) >
                    abs(temp.year - ZXCVBN_DATE_REF_YEAR))
                memcpy(&best, &temp, sizeof(temp));
        }
    }
    if (!best.day)
        return 0;
    memcpy(date, &best, sizeof(best));
    return 1;
}

static uint8_t
zxcvbn_date_probe_split(struct zxcvbn_date *date,
                        char *str, uint32_t len, uint8_t *split,
                        struct zxcvbn_date *dates, unsigned int dates_num)
{
    uint8_t len2, flags;
    uint16_t nums[3];

    len2 = len - split[1];
    nums[0] = zxcvbn_parse_number(str, split[0]);
    nums[1] = zxcvbn_parse_number(str + split[0], split[1] - split[0]);
    nums[2] = zxcvbn_parse_number(str + split[1], len2);

    flags = ZXCVBN_DATE_PROBE_LEFT_YEAR | ZXCVBN_DATE_PROBE_RIGHT_YEAR;
    if (split[0] == 4 || len2 == 1)
        flags &= ~ZXCVBN_DATE_PROBE_RIGHT_YEAR;
    else if (split[0] == 1 || len2 == 4)
        flags &= ~ZXCVBN_DATE_PROBE_LEFT_YEAR;
    if (split[0] == 4 || len2 == 4)
        flags |= ZXCVBN_DATE_PROBE_FULL_YEAR;
    return zxcvbn_date_probe(date, nums, flags, dates, dates_num);
}

static struct zxcvbn_match *
zxcvbn_date_add_match(struct zxcvbn_res *res, uint32_t i, uint32_t j,
                      struct zxcvbn_date *date)
{
    struct zxcvbn_match *match;

    if (!(match = match_add(res)))
        return NULL;

    match->type = ZXCVBN_MATCH_TYPE_DATE;
    match->i = i;
    match->j = j;
    memcpy(&match->date, date, sizeof(*date));
    return match;
}

static int8_t
zxcvbn_date_match_nosep(struct zxcvbn_res *res,
                        char *password, int password_len,
                        struct zxcvbn_date *dates, unsigned int dates_num)
{
    static uint8_t split4[][2] = {{1, 2}, {2, 3}, {0, 0}},
                   split5[][2] = {{1, 3}, {2, 3}, {0, 0}},
                   split6[][2] = {{1, 2}, {2, 4}, {4, 5}, {0, 0}},
                   split7[][2] = {{1, 3}, {2, 3}, {4, 5}, {4, 6}, {0, 0}},
                   split8[][2] = {{2, 4}, {4, 6}, {0, 0}};
    static uint8_t *splits[] = {
        split4[0], split5[0], split6[0], split7[0], split8[0]
    };
    struct zxcvbn_date best, date;
    uint32_t len, i, j, k, d;
    uint8_t *split;

    i = 0;
    while (i + ZXCVBN_DATE_MIN_NOSEP_LEN - 1 < password_len) {
        if (!isdigit(password[i])) {
            i++;
            continue;
        }
        for (j = i + 1; j < password_len && isdigit(password[j]); j++)
            {;}
        len = j - i;
        if (len < ZXCVBN_DATE_MIN_NOSEP_LEN) {
            i += len + 1;
            continue;
        }

        for (j = MIN(len, ZXCVBN_DATE_MAX_NOSEP_LEN);
                                j >= ZXCVBN_DATE_MIN_NOSEP_LEN; j--) {
            for (k = i; k <= i + len - j; k++) {

                /* probe only year */
                if (j == 4 && zxcvbn_date_probe_year(&date, password + k)) {
                    for (d = 0; d < dates_num; d++) {
                        if (date.year == dates[d].year) {
                            date.flags |= ZXCVBN_DATE_FROM_LIST;
                            break;
                        }
                    }
                    if (!zxcvbn_date_add_match(res, k, k + j - 1, &date))
                        return -1;
                    continue;
                }

                /* probe full date */
                memset(&best, 0, sizeof(best));
                split = splits[j - ZXCVBN_DATE_MIN_NOSEP_LEN];
                while (split[0]) {
                    if (zxcvbn_date_probe_split(&date, password + k, j, split,
                                                dates, dates_num)) {
                        if (!best.day ||
                                (date.flags & ZXCVBN_DATE_FROM_LIST) ||
                                abs(best.year - ZXCVBN_DATE_REF_YEAR) >
                                abs(date.year - ZXCVBN_DATE_REF_YEAR))
                            memcpy(&best, &date, sizeof(date));
                        if (best.flags & ZXCVBN_DATE_FROM_LIST)
                            break;
                    }
                    split += 2;
                }
                if (best.day) {
                    if (!zxcvbn_date_add_match(res, k, k + j - 1, &best))
                        return -1;
                }
            }
        }
        i += len + 1;
    }
    return 0;
}

static int8_t
zxcvbn_date_match_sep(struct zxcvbn_res *res,
                      char *password, int password_len,
                      struct zxcvbn_date *dates, unsigned int dates_num)
{
    static struct zxcvbn_date_state states[] = {
        /*          d   s   x       skip  num try p_fl */
        /*  0 */ {{ 1, 15, -1}, { 1,  1,  2}, -1,  0,   0},
        /*  1 */ {{28,  2, -1}, { 1,  1,  3}, -1,  0,   0},
        /*  2 */ {{ 3, -1, -1}, { 1,  4,  4},  0,  0,   0},
        /*  3 */ {{ 4, 10, -1}, { 1,  1,  5}, -1,  0,   0},
        /*  4 */ {{-1,  5, -1}, { 3,  1,  6}, -1,  0,   0},
        /*  5 */ {{ 6, -1, -1}, { 1,  7,  7},  1,  0,   0},
        /*  6 */ {{ 7, -1, -1}, { 1,  3,  8},  2,  1,   1},
        /*  7 */ {{ 8, -1, -1}, { 1,  1,  9},  2,  1,   3},
        /*  8 */ {{ 9, -1, -1}, { 1,  1, 10}, -1,  0,   0},
        /*  9 */ {{-1, -1, -1}, { 1,  1, 11},  2,  1,   6},
        /* 10 */ {{11, -1, -1}, { 1,  6,  6},  1,  0,   0},
        /* 11 */ {{12, -1, -1}, { 1,  3,  7},  2,  1,   1},
        /* 12 */ {{13, -1, -1}, { 1,  1,  8},  2,  1,   3},
        /* 13 */ {{14, -1, -1}, { 1,  1,  9}, -1,  0,   0},
        /* 14 */ {{-1, -1, -1}, { 1,  1, 10},  2,  1,   6},
        /* 15 */ {{16, -1, -1}, { 1,  3,  3},  0,  0,   0},
        /* 16 */ {{17, 23, -1}, { 1,  1,  4}, -1,  0,   0},
        /* 17 */ {{-1, 18, -1}, { 2,  1,  5}, -1,  0,   0},
        /* 18 */ {{19, -1, -1}, { 1,  6,  6},  1,  0,   0},
        /* 19 */ {{20, -1, -1}, { 1,  2,  7}, -1,  0,   0},
        /* 20 */ {{21, -1, -1}, { 1,  2,  8},  2,  1,   2},
        /* 21 */ {{22, -1, -1}, { 1,  6,  9}, -1,  0,   0},
        /* 22 */ {{-1, -1, -1}, { 6,  5, 10},  2,  1,   6},
        /* 23 */ {{24, -1, -1}, { 1,  5,  5},  1,  0,   0},
        /* 24 */ {{25, -1, -1}, { 1,  2,  6}, -1,  0,   0},
        /* 25 */ {{26, -1, -1}, { 1,  2,  7},  2,  1,   2},
        /* 26 */ {{27, -1, -1}, { 1,  5,  8}, -1,  0,   0},
        /* 27 */ {{-1, -1, -1}, { 5,  4,  9},  2,  1,   6},
        /* 28 */ {{29, -1, -1}, { 1,  1,  4}, -1,  0,   0},
        /* 29 */ {{-1, 30, -1}, { 1,  1,  5}, -1,  0,   0},
        /* 30 */ {{31, -1, -1}, { 1,  6,  6},  0,  0,   0},
        /* 31 */ {{32, 36, -1}, { 1,  1,  7}, -1,  0,   0},
        /* 32 */ {{-1, 33, -1}, { 5,  1,  8}, -1,  0,   0},
        /* 33 */ {{34, -1, -1}, { 1,  9,  9},  1,  0,   0},
        /* 34 */ {{35, -1, -1}, { 1,  2, 10},  2,  1,   5},
        /* 35 */ {{-1, -1, -1}, { 2,  2, 11},  2,  1,   5},
        /* 36 */ {{37, -1, -1}, { 1,  8,  8},  1,  0,   0},
        /* 37 */ {{38, -1, -1}, { 1,  2,  9},  2,  1,   5},
        /* 38 */ {{-1, -1, -1}, { 2,  2, 10},  2,  1,   5},
    };
    struct zxcvbn_date best, date;
    struct zxcvbn_date_state *state;
    uint16_t n, nums[3];
    uint32_t i, j, end;
    uint8_t id, skip;
    int8_t next;
    char ch;

    i = 0;
    end = 0;
    while (i + ZXCVBN_DATE_MIN_SEP_LEN - 1 < password_len) {
        if (!isdigit(password[i])) {
            i++;
            continue;
        }

        state = states;
        memset(&best, 0, sizeof(best));
        n = password[i] - '0';
        for (j = i + 1;; j++) {
            if (j < password_len) {
                ch = password[j];
                if (isdigit(ch)) {
                    id = 0;
                    n = n * 10 + ch - '0';
                } else if (strchr("-._/\\", ch))
                    id = 1;
                else
                    id = 2;
            } else
                id = 2;
            next = state->next[id];
            if (next < 0) {
                skip = state->skip[id];
                break;
            }
            state = states + next;
            if (state->num >= 0)
                nums[state->num] = n;
            if (id)
                n = 0;
            if (!state->try)
                continue;
            if (zxcvbn_date_probe(&date, nums, state->probe_flags,
                                  dates, dates_num)) {
                int replace = 0;

                if (!best.day)
                    replace = 1;
                else {
                    if (best.flags & ZXCVBN_DATE_FROM_LIST) {
                        if ((date.flags & ZXCVBN_DATE_FROM_LIST) && end < j)
                            replace = 1;
                    } else {
                        if ((date.flags & ZXCVBN_DATE_FROM_LIST) ||
                                abs(best.year - ZXCVBN_DATE_REF_YEAR) >
                                abs(date.year - ZXCVBN_DATE_REF_YEAR) ||
                                end < j)
                            replace = 1;
                    }
                }
                if (replace) {
                    memcpy(&best, &date, sizeof(date));
                    end = j;
                }
            }
        }
        if (best.day) {
            best.flags |= ZXCVBN_DATE_SEPARATOR;
            if (!zxcvbn_date_add_match(res, i, end, &best))
                return -1;
        }
        i += skip;
    }
    return 0;
}

static int8_t
zxcvbn_date_match(struct zxcvbn_res *res, char *password, int password_len,
                  struct zxcvbn_date *dates, unsigned int dates_num)
{
    if (zxcvbn_date_match_nosep(res, password, password_len, dates, dates_num))
        return -1;
    if (zxcvbn_date_match_sep(res, password, password_len, dates, dates_num))
        return -1;
    return 0;
}

static void
zxcvbn_date_calculate_entropy(struct zxcvbn *zxcvbn, struct zxcvbn_match *match)
{
    double possib;

    if (match->date.flags & ZXCVBN_DATE_FROM_LIST)
        match->entropy = 0;
    else {
        possib = fmax(abs(match->date.year - ZXCVBN_DATE_REF_YEAR),
                             ZXCVBN_DATE_MIN_YEAR_SPACE);
        if (!(match->date.flags & ZXCVBN_DATE_ONLY_YEAR))
            possib *= 12 * 31;
        match->entropy = log2(possib);
    }
    if (match->date.flags & ZXCVBN_DATE_FULL_YEAR)
        match->entropy += 1;
    if (match->date.flags & ZXCVBN_DATE_SEPARATOR)
        match->entropy += 2;
}

/* Date ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

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
        dst[i] = zxcvbn->pack_table[(unsigned char) tolower(src[i])];

    return dst;
}

static int
match_dict(struct zxcvbn_res *res, const char *password, unsigned int password_len,
           char **dict_words, unsigned int n_dict_words)
{
    int i, dict_word_len, remain;
    char pack_password[ZXCVBN_PASSWORD_LEN_MAX], *pack_dict_word, *s;
    struct zxcvbn_dict *dict;
    struct zxcvbn *zxcvbn;

    zxcvbn = res->zxcvbn;

    pack_word(zxcvbn, pack_password, password, password_len);

    for (i = 0; i < n_dict_words; ++i) {
        dict_word_len = strlen(dict_words[i]);
        if (!dict_word_len || password_len < dict_word_len)
            continue;
        pack_dict_word = pack_word(zxcvbn, dict_words[i],
                                   dict_words[i], dict_word_len);
        s = pack_password;
        while (1) {
            if ((remain = password_len - (s - pack_password)) <= 0 ||
                    !(s = memmem(s, remain, pack_dict_word, dict_word_len)))
                break;
            if (!push_match_dict(res, s - pack_password,
                                 s - pack_password + dict_word_len - 1, 1))
                return -1;
            s += dict_word_len;
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
    double possibilities;

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

    if (upper == 1 && isupper(password[match->i]))
        match->entropy += 1;
    else if (upper) {
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
    double possibilities;
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

struct zxcvbn *
zxcvbn_init_ex(struct zxcvbn *zxcvbn, struct zxcvbn_opts *opts)
{
    static struct zxcvbn_opts default_opts;
    zxcvbn_malloc_t malloc_;
    int i, l33t;
    const char *s;

    if (!opts)
        opts = &default_opts;
    malloc_ = (opts->malloc ? opts->malloc : malloc);

    if (zxcvbn)
        memset(zxcvbn, 0, sizeof(*zxcvbn));
    else {
        if (!(zxcvbn = malloc_(sizeof(*zxcvbn))))
            return NULL;
        memset(zxcvbn, 0, sizeof(*zxcvbn));
        zxcvbn->allocated = 1;
    }

    zxcvbn->zxcvbn_malloc = malloc_;
    if (opts->malloc) {
        assert(opts->realloc);
        assert(opts->free);

        zxcvbn->zxcvbn_realloc = opts->realloc;
        zxcvbn->zxcvbn_free = opts->free;
    } else {
        assert(!opts->realloc);
        assert(!opts->free);

        zxcvbn->zxcvbn_realloc = realloc;
        zxcvbn->zxcvbn_free = free;
    }

    zxcvbn->max_matches_num = opts->max_matches_num;

    LIST_INIT(&zxcvbn->dict_head);

    memset(zxcvbn->pack_table, '.', sizeof(zxcvbn->pack_table));

    for (i = 'a'; i <= 'z'; ++i)
        zxcvbn->pack_table[i] = zxcvbn->pack_table_size++;

    for (i = '0'; i <= '9'; ++i)
        zxcvbn->pack_table[i] = zxcvbn->pack_table_size++;

    for (s = opts->symbols; *s != '\0'; ++s) {
        if (zxcvbn->pack_table[*s] == '.') {
            zxcvbn->pack_table[*s] = zxcvbn->pack_table_size++;
            zxcvbn->n_symbols++;
        }
    }

    // l33t
    for (i = 0; i < 256; ++i) { 
        switch (i) {
        case '4':
        case '@':
            l33t = 'a';
            break;
        case '8':
            l33t = 'b';
            break;
        case '(':
        case '{':
        case '[':
        case '<':
            l33t = 'c';
            break;
        case '3':
            l33t = 'e';
            break;
        case '6':
        case '9':
            l33t = 'g';
            break;
        case '1':
        case '!':
        case '|':
            l33t = 'i';
            break;
        case '0':
            l33t = 'o';
            break;
        case '$':
        case '5':
            l33t = 's';
            break;
        case '+':
        case '7':
            l33t = 't';
            break;
        case '%':
            l33t = 'x';
            break;
        case '2':
            l33t = 'z';
            break;
        default:
            continue;
        }
        zxcvbn->pack_table[i] = zxcvbn->pack_table[l33t];
    }
    make_spatial_graph(zxcvbn);

    return zxcvbn;
}

struct zxcvbn *
zxcvbn_init(struct zxcvbn *zxcvbn_buf,
            void *(*zxcvbn_malloc)(size_t size),
            void *(*zxcvbn_realloc)(void *ptr, size_t size),
            void (*zxcvbn_free)(void *ptr),
            const char *symbols)
{
    struct zxcvbn_opts opts;

    zxcvbn_opts_init(&opts);
    opts.malloc = zxcvbn_malloc;
    opts.realloc = zxcvbn_realloc;
    opts.free = zxcvbn_free;
    opts.symbols = symbols;
    return zxcvbn_init_ex(zxcvbn_buf, &opts);
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
zxcvbn_match_ex(struct zxcvbn_res *res,
                const char *password,      unsigned int password_len,
                char **words,              unsigned int words_num,
                struct zxcvbn_date *dates, unsigned int dates_num)
{
    int i;
    struct zxcvbn *zxcvbn;
    struct zxcvbn_match *match;

    assert(password_len > 0);
    assert(password_len <= ZXCVBN_PASSWORD_LEN_MAX);

    if (!dates)
        dates_num = 0;

    if (match_spatial(res, password, password_len))
        return -1;
    if (match_digits(res, password, password_len))
        return -1;
    if (zxcvbn_date_match(res, (char *) password, password_len,
                          dates, dates_num))
        return -1;
    if (zxcvbn_sequence_match(res, (char *) password, password_len))
        return -1;
    if (zxcvbn_repeat_match(res, (char *) password, password_len))
        return -1;
    if (match_dict(res, password, password_len, words, words_num))
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
        case ZXCVBN_MATCH_TYPE_SEQUENCE:
            zxcvbn_sequence_calculate_entropy(zxcvbn, match,
                                              (char *) password, password_len);
            break;
        case ZXCVBN_MATCH_TYPE_REPEAT:
            zxcvbn_repeat_calculate_entropy(zxcvbn, match,
                                            (char *) password, password_len);
            break;
        case ZXCVBN_MATCH_TYPE_DATE:
            zxcvbn_date_calculate_entropy(zxcvbn, match);
            break;
        default:
            assert(0);
        }
    }

    return min_entropy(res, password, password_len);
}

int
zxcvbn_match(struct zxcvbn_res *res,
             const char *password, unsigned int password_len,
             char **words,         unsigned int words_num)
{
    return zxcvbn_match_ex(res, password, password_len, words, words_num,
                           NULL, 0);
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

static void
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
    dict->name[sizeof(dict->name) - 1] = '\0';
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
            if (node->rank == -1 || node->rank > rank)
                node->rank = rank;
            break;
        }

        parent = node;
    }

    return 1;
}
