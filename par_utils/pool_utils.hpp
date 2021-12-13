
int min(int a, int b)
{
  if (a < b)
  {
    return a;
  }
  return b;
}

int max(int a, int b)
{
  if (a > b)
  {
    return a;
  }
  return b;
}


float sum(float *a, int n)
{
  float s = 0.0;
  for (int i = 0; i < n; i++)
  {
    s += a[i];
  }
  return s;
}