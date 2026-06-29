import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.persistence_air_coral.ismPersistenceAirCoral;
import model.ModelExecutionTest;

public class PersistenceAirCoralTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismPersistenceAirCoral();
    }
}
