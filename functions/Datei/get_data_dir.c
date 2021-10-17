#include "config.h"
#include "utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sys/stat.h>

#define CARAT_PATH_MAX 4096

// location of the current executable
static char carat_exec_location[CARAT_PATH_MAX] = "";
static char carat_data_location[CARAT_PATH_MAX] = "";

// The function 'find_yourself' is based on code from GAP, which in turn
// is based on code (C) 2015 Mark Whitis, under the MIT License :
// https://stackoverflow.com/a/34271901/928031
static void
find_yourself(const char * argv0, char * result, size_t resultsize)
{
    assert(resultsize >= CARAT_PATH_MAX);

    char tmpbuf[CARAT_PATH_MAX];

    // absolute path, like '/usr/bin/gap'
    if (argv0[0] == '/') {
        if (realpath(argv0, result) && !access(result, F_OK)) {
            return;    // success
        }
    }
    // relative path, like 'bin/gap.sh'
    else if (strchr(argv0, '/')) {
        if (!getcwd(tmpbuf, sizeof(tmpbuf)))
            return;
        strlcat(tmpbuf, "/", sizeof(tmpbuf));
        strlcat(tmpbuf, argv0, sizeof(tmpbuf));
        if (realpath(tmpbuf, result) && !access(result, F_OK)) {
            return;    // success
        }
    }
    // executable name, like 'gap'
    else {
        char pathenv[CARAT_PATH_MAX], *saveptr, *pathitem;
        strlcpy(pathenv, getenv("PATH"), sizeof(pathenv));
        pathitem = strtok_r(pathenv, ":", &saveptr);
        for (; pathitem; pathitem = strtok_r(NULL, ":", &saveptr)) {
            strlcpy(tmpbuf, pathitem, sizeof(tmpbuf));
            strlcat(tmpbuf, "/", sizeof(tmpbuf));
            strlcat(tmpbuf, argv0, sizeof(tmpbuf));
            if (realpath(tmpbuf, result) && !access(result, F_OK)) {
                return;    // success
            }
        }
    }

    *result = 0;    // reset buffer after error
}

// the following is directly based on GAP's SetupGAPLocation
void setup_carat_location(const char * argv0)
{
    // In the code below, we keep reseting locBuf, as some of the methods we
    // try do not promise to leave the buffer empty on a failed return.
    char locBuf[CARAT_PATH_MAX] = "";

#if defined(__APPLE__) && defined(__MACH__)
    uint32_t len = sizeof(locBuf);
    if (_NSGetExecutablePath(locBuf, &len) != 0) {
        *locBuf = 0;    // reset buffer after error
    }
#endif

    // try Linux procfs
    if (!*locBuf) {
        ssize_t ret = readlink("/proc/self/exe", locBuf, sizeof(locBuf));
        if (ret < 0)
            *locBuf = 0;    // reset buffer after error
    }

    // try FreeBSD / DragonFly BSD procfs
    if (!*locBuf) {
        ssize_t ret = readlink("/proc/curproc/file", locBuf, sizeof(locBuf));
        if (ret < 0)
            *locBuf = 0;    // reset buffer after error
    }

    // try NetBSD procfs
    if (!*locBuf) {
        ssize_t ret = readlink("/proc/curproc/exe", locBuf, sizeof(locBuf));
        if (ret < 0)
            *locBuf = 0;    // reset buffer after error
    }

    // if we are still failing, go and search the path
    if (!*locBuf) {
        find_yourself(argv0, locBuf, CARAT_PATH_MAX);
    }

    // resolve symlinks (if present)
    if (!realpath(locBuf, carat_exec_location))
        *carat_exec_location = 0;    // reset buffer after error

    // now strip the executable name off
    int length = strlen(carat_exec_location);
    while (length > 0 && carat_exec_location[length] != '/') {
        carat_exec_location[length] = 0;
        length--;
    }
}

static int is_valid_dir(const char * path)
{
    struct stat buf;
    if (!path || !*path)
        return 0;
    if (stat(path, &buf))
        return 0;
    return S_ISDIR(buf.st_mode);
}

static void setup_data_dir(void)
{
    char * dir = getenv("CARAT_DIR");
    if (dir) {
        strcpy(carat_data_location, dir);
        strcat(carat_data_location, "/tables");
        if (is_valid_dir(carat_data_location)) {
            return;
        }
    }

    // perhaps this is a binary which was never moved from its compile location?
    strcpy(carat_data_location, carat_exec_location);
    strcat(carat_data_location, "/../tables");
    if (is_valid_dir(carat_data_location)) {
        return;
    }

    // perhaps this is a binary which was installed into PREFIX/bin, and the
    // data tables were installed into RPEFIX/share/carat/tables?
    strcpy(carat_data_location, carat_exec_location);
    strcat(carat_data_location, "/../share/carat/tables");
    if (is_valid_dir(carat_data_location)) {
        return;
    }

    // instead of looking relative to the executable, try
    // RPEFIX/share/carat/tables
    strcpy(carat_data_location, DATADIR "/carat/tables");
    if (is_valid_dir(carat_data_location)) {
        return;
    }

    // giving up
    fprintf(stderr, "could not locate carat data directory");
    exit(1);
}

const char * get_data_dir(void)
{
    if (carat_data_location[0] == '\0')
        setup_data_dir();
    return carat_data_location;
}
