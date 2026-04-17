/* Miscellaneous utilities for workflows from the ONT Customer Workflows Group.
 */
class CWUtil {

    /* Mutate the global Nextflow params map
    *
    * Occasionally, we may wish to mutate the value of a parameter provided
    * by the user. Typically, this leads to workflows with `params.my_param`
    * and `params._my_param` which is ripe for confusion. Instead, we can
    * mutate the parameter value in the Nextflow params ScriptMap itself
    * with the following call:
    *
    *     CWUtil.mutateParam(params, k, v)
    *
    * This is possible as Groovy actually has a surprisingly loose
    * definition of "private", and allows us to call the private `allowNames`
    * method on the ScriptMap which removes the read-only status for a key set.
    * We can follow this up with a call to the private `put0` to reinsert
    * the key and mark it as read-only again.
    */
    public static void mutateParam(nf_params, key, value) {
        Set s = [key] // must be a set to allow call to allowNames
        nf_params.allowNames(s)
        nf_params.put0(key, value)
    }

    /* Safely get the workflow entrypoint from the Nextflow params map 
    *
    * Our workflows support alternative "entrypoints". For example, a workflow
    * may expose an entrypoint that can be used to validate parameters and
    * dry run the workflow, or to automatically set different defaults or
    * behaviours. This is controlled by setting `params.wf.entrypoint`.
    *
    * This util safely checks and extracts this parameter if set and avoids
    * "Access to undefined parameter" errors otherwise raised during strict
    * Nextflow execution.
    * 
    * Note that entrypoint is expected to be null for the "main" entrypoint.
    */
    public static def getEntrypoint(nf_params) {
        def entrypoint = null
        if( nf_params.containsKey('wf')
            && nf_params.wf instanceof Map
            && nf_params.wf.containsKey('entrypoint') )
        {
            entrypoint = nf_params.wf.entrypoint
        }
        return entrypoint
    }
}
