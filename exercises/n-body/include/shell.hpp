#ifndef SHELL_H_
#define SHELL_H_

class Shell {
  public:
    Shell();
    Shell(Shell &&) = default;
    Shell(const Shell &) = default;
    Shell &operator=(Shell &&) = default;
    Shell &operator=(const Shell &) = default;
    ~Shell();

  private:
};

#endif // ! SHELL_H_
