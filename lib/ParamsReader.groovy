//package utils

import groovy.json.JsonSlurper

class ParamsReader {

    public static Object readParamsFromJsonSettings(String path) throws Exception {
        return new JsonSlurper().parseText(new File(path).text)
    }

}