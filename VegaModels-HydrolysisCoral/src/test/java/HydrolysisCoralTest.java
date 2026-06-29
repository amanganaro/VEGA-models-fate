import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.hydrolysis_coral.ismHydrolysisCoral;
import model.ModelExecutionTest;

public class HydrolysisCoralTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismHydrolysisCoral();
    }
}
