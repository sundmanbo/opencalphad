/*
 *      @(#) Driver for reading a character from keyboard in raw I/O mode
 */
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <sys/ioctl.h>

/*--------------------------*/
#ifdef BSD
#include <sgtty.h>
#else
/*--------------------------*/
#ifdef Linux
#define G77
#endif

#ifdef CYGWIN
#define G77
#endif

#ifdef G77
#include <termio.h>
#else
#include <sys/termio.h>
#endif
/*--------------------------*/
#endif

#include <signal.h>

/*V13 #include <sys/types.h> */
/*V13 #include <termios.h> */
/******************************************************************************/
/* return the next key typed in hot (raw I/O) mode.  */
char getkeyC(void) {
#ifdef BSD
        struct sgttyb   oldtty, newtty;
        char            c;

        ioctl(0, TIOCGETP, &oldtty);

        newtty = oldtty;
        newtty.sg_flags = RAW;

        ioctl(0, TIOCSETP, &newtty);

        read(0, &c, 1);

        ioctl(0, TIOCSETP, &oldtty);
#else
        struct termio   oldtty, newtty;
/*V13   struct termios   oldtty, newtty; */
        char            c;

        ioctl(0, TCGETA, &oldtty);
/*V13 	tcgetattr(0, &oldtty);  */

        newtty = oldtty;
        newtty.c_iflag = BRKINT | IXON | ISTRIP;
        newtty.c_lflag = 0;
        newtty.c_cc[VEOF] = 1;

        ioctl(0, TCSETA, &newtty);
/*V13 	tcsetattr(0, TCSANOW, &newtty); */

        read(0, &c, 1);

        ioctl(0, TCSETA, &oldtty);
/*V13 	tcsetattr(0, TCSANOW, &oldtty); */
#endif
        /* fprintf(stderr,"C:c=%c\n",c); */
	/* fflush(stdout); */
        return(c);
}
/******************************************************************************/
/* Commonly, a C routine called name_ can be called from Fortran as name; plus less-common ones */
int getkey4f_(void) { return(getkeyC()); }
int _getkey4f(void) { return(getkeyC()); }
int getkey4f(void) { return(getkeyC()); }
int GETKEY4F(void)  { return(getkeyC()); }


/* http://www.urbanjost.altervista.org/LIBRARY/libCLI/Getkey/getkey.html */
