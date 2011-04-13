
class Options {
public:
    Options()
        :  file_input(false), draw(false), statistics(false),
           check(false), number_of_points(0),
           min(-1.1), max(1.1), winx(400), winy(400)
    {}

    char program[100];
    char fname[100];
    bool file_input;
    bool draw;
    bool statistics;
    bool check;
    int  number_of_points;
    double min;
    double max;
    int winx;
    int winy;
};


void usage(char* program);

bool parse(int argc, char* argv[], Options &opt);
