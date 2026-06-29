import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.persistence_quantative_water_irfmn.ismPersistenceWaterQuantitativeIrfmn;
import model.ModelExecutionTest;

public class PersistenceWaterQuantitativeIRFMNTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismPersistenceWaterQuantitativeIrfmn();
    }
}
