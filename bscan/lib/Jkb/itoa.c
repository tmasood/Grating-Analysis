void reverse(char *);

void itoa(int n, char *s)
{
  int i, sign;

  if ((sign = n) < 0) /* record sign */
    n = -n;
  i = 0;
  do {     /* generate in reverse order */
    s[i++] = n % 10 + '0'; /* get next digit */
    } while ((n /= 10) > 0); /* delete it */
  if ( sign < 0 )
    s[i++] = '-';
  s[i] = '\0';
  reverse(s);
}
