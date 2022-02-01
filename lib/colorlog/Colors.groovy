package colorlog

import groovy.transform.CompileStatic


@CompileStatic
interface Colors {
    // Reset
    final static RESET = "\033[0m";

    // Dim
    final static String DIM = "\033[2m";

    // Regular Colors
    final static String BLACK = "\033[0;30m";
    final static String RED = "\033[0;31m";
    final static String GREEN = "\033[0;32m";
    final static String YELLOW = "\033[0;33m";
    final static String BLUE = "\033[0;34m";
    final static String PURPLE = "\033[0;35m";
    final static String CYAN = "\033[0;36m";
    final static String WHITE = "\033[0;37m";

    // Bold
    final static String BLACK_BOLD = "\033[1;30m";
    final static String RED_BOLD = "\033[1;31m";
    final static String GREEN_BOLD = "\033[1;32m";
    final static String YELLOW_BOLD = "\033[1;33m";
    final static String BLUE_BOLD = "\033[1;34m";
    final static String PURPLE_BOLD = "\033[1;35m";
    final static String CYAN_BOLD = "\033[1;36m";
    final static String WHITE_BOLD = "\033[1;37m";

    // Underline
    final static String BLACK_UNDERLINED = "\033[4;30m";
    final static String RED_UNDERLINED = "\033[4;31m";
    final static String GREEN_UNDERLINED = "\033[4;32m";
    final static String YELLOW_UNDERLINED = "\033[4;33m";
    final static String BLUE_UNDERLINED = "\033[4;34m";
    final static String PURPLE_UNDERLINED = "\033[4;35m";
    final static String CYAN_UNDERLINED = "\033[4;36m";
    final static String WHITE_UNDERLINED = "\033[4;37m";

    // Background
    final static String BLACK_BACKGROUND = "\033[40m";
    final static String RED_BACKGROUND = "\033[41m";
    final static String GREEN_BACKGROUND = "\033[42m";
    final static String YELLOW_BACKGROUND = "\033[43m";
    final static String BLUE_BACKGROUND = "\033[44m";
    final static String PURPLE_BACKGROUND = "\033[45m";
    final static String CYAN_BACKGROUND = "\033[46m";
    final static String WHITE_BACKGROUND = "\033[47m";

    // High Intensity
    final static  String BLACK_BRIGHT = "\033[0;90m";
    final static  String RED_BRIGHT = "\033[0;91m";
    final static  String GREEN_BRIGHT = "\033[0;92m";
    final static  String YELLOW_BRIGHT = "\033[0;93m";
    final static  String BLUE_BRIGHT = "\033[0;94m";
    final static  String PURPLE_BRIGHT = "\033[0;95m";
    final static  String CYAN_BRIGHT = "\033[0;96m";
    final static  String WHITE_BRIGHT = "\033[0;97m";

    // Bold High Intensity
    final static  String BLACK_BOLD_BRIGHT = "\033[1;90m";
    final static  String RED_BOLD_BRIGHT = "\033[1;91m";
    final static  String GREEN_BOLD_BRIGHT = "\033[1;92m";
    final static  String YELLOW_BOLD_BRIGHT = "\033[1;93m";
    final static  String BLUE_BOLD_BRIGHT = "\033[1;94m";
    final static  String PURPLE_BOLD_BRIGHT = "\033[1;95m";
    final static  String CYAN_BOLD_BRIGHT = "\033[1;96m";
    final static  String WHITE_BOLD_BRIGHT = "\033[1;97m";

    // High Intensity backgrounds
    final static  String BLACK_BACKGROUND_BRIGHT = "\033[0;100m";
    final static  String RED_BACKGROUND_BRIGHT = "\033[0;101m";
    final static  String GREEN_BACKGROUND_BRIGHT = "\033[0;102m";
    final static  String YELLOW_BACKGROUND_BRIGHT = "\033[0;103m";
    final static  String BLUE_BACKGROUND_BRIGHT = "\033[0;104m";
    final static  String PURPLE_BACKGROUND_BRIGHT = "\033[0;105m";
    final static  String CYAN_BACKGROUND_BRIGHT = "\033[0;106m";
    final static  String WHITE_BACKGROUND_BRIGHT = "\033[0;107m";
}