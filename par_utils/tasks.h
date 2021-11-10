class Task
{
public:
  virtual void execute() {};
};

class Dot_task : public Task
{
private:
  int start;
  int stop;
  float *v1;
  float *v2;
  float *res_ptr;

public:
  Dot_task(int s, int st, float *x, float *y, float *r);
  void execute();
};
