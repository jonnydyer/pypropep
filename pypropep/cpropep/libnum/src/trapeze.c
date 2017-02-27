
int trapeze(double *data, int n_point, int col, int off, double *integral)
{
  int i;
  double val = 0.0;

  for (i = 0; i < n_point - 1; i++)
  {
    val += (data[off + i*col] + data[off + (i + 1)*col])*
      (data[0 + (i + 1)*col] - data[0 + i*col])/2;
  }

  *integral = val;

  return 0;
}
