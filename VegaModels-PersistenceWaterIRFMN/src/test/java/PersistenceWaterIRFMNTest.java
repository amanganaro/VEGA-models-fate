import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.persistence_water_irfmn.ismPersistenceWaterIrfmn;
import model.ModelExecutionTest;

public class PersistenceWaterIRFMNTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismPersistenceWaterIrfmn();
    }
}
