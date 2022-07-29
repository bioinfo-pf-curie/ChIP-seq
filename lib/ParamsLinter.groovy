//package utils

import org.slf4j.Logger

// TODO: Check updates for the new nf core JSON schema
// https://github.com/nf-core/tools/issues/687
class ParamsLinter {

    LinkedHashMap params
    Object paramsWithUsage
    public static Logger log

    ParamsLinter(params, paramsWithUsage, log){
        this.params = params
        this.paramsWithUsage = paramsWithUsage
        this.log = log

        check()
    }

    LinkedHashMap lint() {
        // Lower case and split params with choices
        def newParams = this.params.each {param ->
            def usage = this.paramsWithUsage.find{us -> us.name == param.key}
            if (param.value && usage && usage.choices ) {
                this.params[param.key] = param.value.toString().contains(',') ? formatParamList(param.value) : param.value == "all" ? usage.choices: param.value
            }
        } as LinkedHashMap
        checkUsage(newParams, this.paramsWithUsage)
        return newParams
    }

    void check() {
        checkGenomes(this.params)
        checkUsage(this.params, this.paramsWithUsage)
    }

    static void checkGenomes(params) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
            System.exit(1)
        }
    }

    static void checkUsage(params, paramsWithUsage) {
        paramsWithUsage.each { usage ->
            def paramValue = params.get(usage.name)
            if (paramValue && usage.choices && paramValue instanceof String) {
                if (paramValue.contains(',')){
                    // Exit if option accept only one value
                    if (usage.nargs == 1) {
                        log.error "You can choose only one ${usage.name}, see --help for more information"
                        System.exit(1)
                    }
                    checkParameterList(formatParamList(paramValue), usage.choices)
                } else {
                    checkParameterExistence(paramValue, usage.choices)
                }
            }
        }
    }

    static void checkParameterExistence(it, list) {
        if (!list.contains(it)) {
            log.error "Unknown parameter: ${it} does not belong to the allowed list ${list.toString()}"
            System.exit(1)
        }
    }

    // Compare each parameter with a list of parameters
    static void checkParameterList(list, realList) {
        list.each{
            checkParameterExistence(it, realList)
        }
    }

    static formatParam(param) {
        return param.trim()
    }

    static formatParamList(param) {
        return param.split(',').collect{formatParam(it)}
    }
}
