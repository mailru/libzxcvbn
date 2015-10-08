#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <sys/time.h>
#include "zxcvbn.h"

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(a) (sizeof(a) / sizeof((a)[0]))
#endif

static char *
trim(char *string)
{
    char *b, *e;

    b = string;

    while (*b != '\0') {
        if (strchr(" \r\n\t", *b) == NULL)
            break;
        ++b;
    }

    e = b;

    while (*e != '\0') {
        if (strchr(" \r\n\t", *e) != NULL)
            break;
        ++e;
    }

    string[e - b] = '\0';

    return b;
}

static struct zxcvbn_dict *
read_ranked(struct zxcvbn *zxcvbn, struct zxcvbn_dict *dict_buf, const char *name, const char *path)
{
    unsigned int rank, word_len;
    char word_buf[1024], *word;
    FILE *file;
    struct zxcvbn_dict *dict;

    if ((dict = zxcvbn_dict_init(zxcvbn, dict_buf, name)) == NULL) {
        fprintf(stderr, "zxcvbn_dict_init(\"%s\") failed\n", name);
        return NULL;
    }

    if ((file = fopen(path, "r")) == NULL) {
        fprintf(stderr, "fopen(\"%s\") failed (%d:%s)\n", path, errno, strerror(errno));
        return NULL;
    }

    rank = 1;

    while (fgets(word_buf, sizeof(word_buf), file) != NULL) {
        word = trim(word_buf);
        if ((word_len = strlen(word)) > 0) {
            if (zxcvbn_dict_add_word(dict, word, word_len, rank) < 0) {
                goto err;
            }
                
            ++rank;
        }
    }

    fclose(file);
    return dict;

err:
    fprintf(stderr, "zxcvbn_dict_add_word(\"%s\", \"%.*s\") failed\n", name, word_len, word);
    fclose(file);
    return NULL;
}

static void
print_usage()
{
    printf("Usage: zxcvbn_cli [ -h ] [ -d \"word0 word1 ... wordN\" ] { password0 } [ password1 ] ... [ passwordN]\n");
}

static char *
escape_quotes(char *str)
{
    static char buf[128];
    unsigned int i;

    for (i = 0; *str && i < sizeof(buf) - 1; str++) {
        if (strchr("\"\\", *str))
            buf[i++] = '\\';
        buf[i++] = *str;
    }
    buf[i] = '\0';
    return buf;
}

static void
process_bulk(int argc, char **argv)
{
    char buf[1024], *words[256], *p;
    unsigned int words_num;
    struct zxcvbn_res res;
    struct zxcvbn *z;
    size_t len;
    int opt;

    if (!(z = zxcvbn_init(NULL, NULL, NULL, NULL,
                          "!@#$%^&*()-_+=;:,./?\\|`~[]{}"))) {
        fprintf(stderr, "zxcvbn_init() failed\n");
        exit(EXIT_FAILURE);
    }

    optind = 1;
    opterr = 0;
    while ((opt = getopt(argc, argv, "D:")) != -1) {
        switch (opt) {
            case 'D':
                if (!read_ranked(z, NULL, optarg, optarg))
                    exit(EXIT_FAILURE);
                break;
        }
    }

    while (fgets(buf, sizeof(buf), stdin)) {
        words_num = 0;
        len = strlen(buf);
        if (buf[len - 1] == '\n')
            buf[len - 1] = '\0';
        p = strchr(buf, ' ');
        if (p) {
            *p = '\0';
            for (p = strtok(p + 1, " "); p; p = strtok(NULL, " ")) {
                if (words_num == ARRAY_SIZE(words))
                    break;
                words[words_num++] = p;
            }
        }

        zxcvbn_res_init(&res, z);
        if (zxcvbn_match(&res, buf, strlen(buf),
                         words, words_num) < 0) {
            printf("{\"password\": \"%s\", \"error\": true}\n",
                   escape_quotes(buf));
            fprintf(stderr, "zxcvbn_match(\"%s\") failed\n",
                    escape_quotes(buf));
            zxcvbn_res_release(&res);
            continue;
        }
        printf("{\"password\": \"%s\", \"entropy\": %.1lf}\n",
               escape_quotes(buf), res.entropy);
        zxcvbn_res_release(&res);
    }
    if (ferror(stdin)) {
        fprintf(stderr, "fgets(stdin) failed\n");
        exit(EXIT_FAILURE);
    }
    exit(EXIT_SUCCESS);
}

int
main(int argc, char **argv)
{
    const char *password;
    char *dict_words[256], *dict_word;
    unsigned int n_dict_words;
    int i, opt;
    struct timeval tv0, tv1;
    struct zxcvbn zxcvbn_buf;
    struct zxcvbn *zxcvbn;
    struct zxcvbn_res res;
    struct zxcvbn_match *match;

    n_dict_words = 0;

    if ((zxcvbn = zxcvbn_init(&zxcvbn_buf, NULL, NULL, NULL, "!@#$%^&*()-_+=;:,./?\\|`~[]{}")) == NULL) {
        fprintf(stderr, "zxcvbn_init() failed\n");
        return EXIT_FAILURE;
    }

    while ((opt = getopt(argc, argv, "D:hd:b")) != -1) {
        switch (opt) {
        case 'h':
            print_usage();
            return EXIT_SUCCESS;

        case 'd':
            for (dict_word = strtok(optarg, " "); dict_word != NULL; dict_word = strtok(NULL, " ")) {
                if (n_dict_words == ARRAY_SIZE(dict_words))
                    break;
                dict_words[n_dict_words++] = dict_word;
            }
            break;

        case 'b':
            process_bulk(argc, argv);
            return EXIT_SUCCESS;

        case 'D':
            read_ranked(zxcvbn, NULL, optarg, optarg);
            break;

        default:
            print_usage();
            return EXIT_FAILURE;
        }
    }


    if (optind == argc) {
        print_usage();
        return EXIT_FAILURE;
    }

    for (i = optind; argv[i] != NULL; ++i) {
        password = argv[i];

        zxcvbn_res_init(&res, zxcvbn);

        gettimeofday(&tv0, NULL);
        if (zxcvbn_match(&res, password, strlen(password), dict_words, n_dict_words) < 0) {
            fprintf(stderr, "zxcvbn_match(\"%s\") failed\n", password);
            continue;
        }
        
        gettimeofday(&tv1, NULL);

        printf("t:%lu us\n", (tv1.tv_sec - tv0.tv_sec) * 1000000 + (tv1.tv_usec - tv0.tv_usec));

        printf("password: %s\n", password);
        printf("entropy: %f\n", res.entropy);
        CIRCLEQ_FOREACH(match, &res.match_head, list) {
            printf("\t%s: %.*s -- %f\n", zxcvbn_match_type_string(match->type),
                   match->j - match->i + 1, password + match->i, match->entropy);
        }
        printf("\n");

        zxcvbn_res_release(&res);
    }

    zxcvbn_release(zxcvbn);

    return EXIT_SUCCESS;
}
